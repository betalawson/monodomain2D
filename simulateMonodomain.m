function simulateMonodomain(problem)
% This function simulates the monodomain equation on a mesh where the
% volume fraction (fraction of accessible material, as opposed to say
% collagenous occlusions due to fibrosis) and conductivity tensor are
% allowed to vary between different elements. Completely non-conductive
% elements may also be included in the mesh.
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
stim_times2 = [720];                    % Vector of times to stimulate sites marked as a secondary stimulus (ms)

% Timestepping and solution methods
t_end = 1000;                             % Simulation time (ms)
dt = 0.05;                                % Timestep (ms)
solve_exact = 0;                          % If true, requires exact solves (direct methods) for the linear system solves involved in taking timesteps
preconditioning = 1;                      % If solve_exact = 0, specifies whether to use basic preconditioning (ILU(0)) or not
second_order = 0;                         % Uses second order timestepping. Threatens stability, but provides better accuracy for sufficiently low timestep
lumping_factor = 0;                       % Specifies the amount of mass-lumping to use, in which the mass matrix is relaxed towards a diagonal matrix
                                          % 0 - no lumping, mass matrix comes from integration of a linear interpolant over control volume
                                          % 1 - full lumping, mass matrix is diagonal, by taking row sums of each row of the original mass matrix

% Plotting
visualise = 1;                            % Flag for whether to visualise or not
save_anim = 0;                            % Flag for whether or not to save an animation (filename same as problem name, CAREFUL not to overwite!)
plot_interval = 2;                        % Time interval for plotting (ms)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the scale factor for diffusion (in monodomain model), denoted
% alpha to match paper
alpha = lambda / (lambda + 1) / chi / Cm;

% Perform the encoding of the problem (stiffness and mass matrices) that
% only need to be calculated once. This code also outputs some extra
% calculated mesh information, including a list of which nodes are active
[K, M, mesh] = encodeProblem(problem.occ_map, problem.D_tensor, problem.Vfrac, problem.grid, alpha);

% Finish preparing the numerical method in terms of these matrices
[A_new, A_old, A_J] = prepareNumerics(K, M, dt, second_order, lumping_factor);

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

% If no extra parameters were provided in the problem structure, give them
% a dummy value
if isfield(problem,'extra_params')
    extra_params = problem.extra_params;
else
    extra_params = [];
end

% Create video object if one is needed
if save_anim && visualise
    vid_obj = VideoWriter('outputVideo.avi');
    open(vid_obj);
end

% Initialise problem
[V, S] = initialiseProblem(cell_models, model_assignments, active);

% The old value for the state variables is also set to the current value
% for the first step. Also, the information that comes from the cell model
% (gating variable rate constants and steady states) is initialised as
% blank to show there is no old information for these
b_old = [];
Sinf = [];
invtau = [];
I_stim_old = zeros(N,1);
J_old = [];


% Loop over time integrations
t = 0;
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
    [I_ion, S_new, Sinf, invtau, b_old] = processReaction(V, S, Sinf, invtau, b_old, dt, I_stim, I_stim_old, cell_models, model_assignments, mesh, second_order, extra_params);
    
    % Calculate the total 'current density'
    J = (1/Cm) * ( I_ion(active) + I_stim(active) );
        
    %%% PROCESS UPDATE
    V_active = takeTimestep(V(active), J, J_old, A_new, A_old, A_J, solve_exact, preconditioning, second_order);
    V(active) = V_active;
    
    % Now update the stored current and old values for different variables
    S = S_new;
    I_stim_old = I_stim;
    J_old = J;
    
    % Check plot frequency, plot if hit
    if visualise && ( t - floor(t / plot_interval) * plot_interval <= dt )
        
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
    
end

% Close video object if one was created
if save_anim && visualise
    close(vid_obj);
end