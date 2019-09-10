function CV = calculateCV(problem)
% This function simulates the monodomain equation on a mesh where the
% volume fraction (fraction of accessible material, as opposed to say
% collagenous occlusions due to fibrosis) and conductivity tensor are
% allowed to vary between different elements. Completely non-conductive
% elements may also be included in the mesh. The function measures
% wavespeed, assuming a wave moving horizontally across the domain from the
% left.
%
% A vertex-centred finite volume method is used to numerically integrate
% the monodomain equation, with bilinear interpolation allowing for
% integration of non-diagonal conductivity tensors, as well as non-local
% integration of the other terms.
%
% Input is the problem structure, as created by one of the createProblem
% type files.

% Monodomain parameters
lambda = 1;                               % Ratio of intracellular and extracellular conducitivies
chi = 2000;                               % Surface-to-volume ratio for tissue (cm^-1)
Cm = 1;                                   % Tissue capacitance per unit area (uF/cm²) - cell capacitance is defined in ionic model files

% Stimulus settings
stim_dur = 1;                             % Stimulus duration (ms)
stim_amp = 52;                            % Amplitude of stimulus per unit area (uA/cm²)
stim_times1 = [20];                       % Vector of times to stimulate sites marked as a primary stimulus (ms)
stim_times2 = [352.5];                    % Vector of times to stimulate sites marked as a secondary stimulus (ms)

% Timestepping and solution methods
t_end = 2000;                             % Simulation time (ms)
dt = 0.1;                                % Timestep (ms)
reac_per_diffuse = 1;                     % Number of reaction steps to perform before each diffusive update
timestepping = 'implicit';                % Timestepping for diffusive updates: set to 'implicit' or 'crank-nicholson'
%nonlocal_timederiv = 1;                   % Specifies whether or not to nonlocally integrate the time derivative term
%nonlocal_react = 1;                       % Specifies whether or not to nonlocally integrate the reaction term
solve_exact = 0;                          % Require exact solves (direct methods) for the linear systems that result from the time and space discretisations
react_parallel = 0;                       % Make use of multiple cores in processing of the reaction terms

% Plotting
visualise = 0;                            % Flag for whether to visualise or not
save_anim = 1;                            % Flag for whether or not to save an animation (filename same as problem name, CAREFUL not to overwite!)
plot_interval = 2.5;                      % Time interval for plotting (ms)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the scale factor for diffusion (in monodomain model), denoted
% alpha to match paper
alpha = lambda / (lambda + 1) / chi / Cm;

% Perform the encoding of the problem (stiffness and mass matrices) that
% only need to be calculated once. This code also outputs some extra
% calculated mesh information, including a list of which nodes are active
[K, M, mesh] = encodeProblem(problem.occ_map, problem.D_tensor, problem.Vfrac, problem.grid, alpha);

% Finish preparing the numerical method in terms of these matrices
[A_new, A_old, A_J] = prepareNumerics(K, M, dt, timestepping, nonlocal_timederiv, nonlocal_react);

% Read out the list of active nodes from mesh file for notational
% cleanliness
active = mesh.active;

% Read out the number of nodes to solve at
N = length(active);

% Read out stimulus sites, and convert to a vectorised form
stim_sites1 = problem.stim_sites1;
stim_sites1 = stim_sites1';
stim_sites1 = stim_sites1(:);
stim_sites2 = problem.stim_sites2;
stim_sites2 = stim_sites2';
stim_sites2 = stim_sites2(:);

% Read out cell models, and convert to a vectorised form
cell_models = problem.cell_models;
model_assignments = problem.model_assignments;
model_assignments = model_assignments';
model_assignments = model_assignments(:);

% Read out node X and Y locations, and also vectorise them
nodeX = problem.nodeX;
nodeY = problem.nodeY;
nodeX = nodeX'; nodeX = nodeX(:);
nodeY = nodeY'; nodeY = nodeY(:);

% Find mesh dimensions
Lx = max(nodeX(:)); Ly = max(nodeY(:));

% Set up the wavespeed measurement regions
start_sites = find( nodeX(:) >= 0.35 * Lx );
end_sites = find( nodeX(:) >= 0.95 * Lx );
measure_dist = min( nodeX(end_sites) ) - min( nodeX(start_sites) );
stim_times2 = []; % Turn off other stimulus to be safe


% Create video object if one is needed
if save_anim && visualise
    vid_obj = VideoWriter('outputVideo.avi');
    open(vid_obj);
end

% Initialise problem
[V, S] = initialiseProblem(cell_models, model_assignments, active);

% Loop over time integrations
t = 0;
start_reached = 0; end_reached = 0;
while t < t_end
    
    % Increment time
    t = t + dt;
    
    
    %%% REACTION STEP HANDLING
    
    % Stimulate if this is a stimulus time
    I_stim = zeros(N,1);
    if any( (t - stim_times1) <= stim_dur & (t - stim_times1) >= 0 )
        I_stim(stim_sites1) = -stim_amp;
    end
    
    if any( (t - stim_times2) <= stim_dur & (t - stim_times2) >= 0 )
        I_stim(stim_sites2) = -stim_amp;
    end
    
    % Process reaction update - uses current voltage values and current
    % state variable values, S. Only processes active sites
    [I_ion, S] = processReaction(V, S, dt, I_stim, cell_models, model_assignments, mesh, react_parallel);
    
    % Calculate the total 'current density'
    J = (1/Cm) * ( I_ion(active) + I_stim(active) );
    
    
    %%% PROCESS UPDATE
    
    V_active = takeTimestep(V(active), J, A_new, A_old, A_J, solve_exact);
    V(active) = V_active;
    
    
    %%% VISUALISATION
    
    % Check plot frequency, plot if hit
    if visualise && ( t - floor(t / plot_interval) * plot_interval <= dt*reac_per_diffuse )
        
        % Visualise the current state
        visualiseState( V, problem.Nx, problem.Ny, problem.occ_map, t );
        
        % Write frame to video object if animation requested
        if save_anim
            frame = getframe(gcf);
            writeVideo(vid_obj, frame);
            % Otherwise, draw on screen so user can watch in real time
        else
            drawnow;
        end
        
    end
    
    
    %%% CONDUCTION VELOCITY CALCULATION
    
    % Check if any of the start sites has been stimulated
    if ~start_reached && any(V(start_sites) > -30)
        start_reached = 1;
        start_time = t;
    end
    
    % Check if the end site has been stimulated
    if ~end_reached && any(V(end_sites) > -30)
        end_reached = 1;
        CV = measure_dist / ((t-start_time) / 1000);
        fprintf('Conduction velocity was %g cm/s \n\n', CV);
        return;
    end
    
    % Check if no sites are active post-stimulus, indicating a wave that
    % died out
    if all(V < -65) && t > 100
        fprintf('Wavefront died out before CV was measured! Fibrosis must have blocked its progress. \n');
        CV = NaN;
        return;
    end
    
end

% Close video object if one was created
if save_anim && visualise
    close(vid_obj);
end

% Program should exit after CV was successfully calculated. Output a
% message here indicating something else going wrong
fprintf('CV measurement failed. Excitation still present, implying a longer simulation time was required. \n');
CV = NaN;