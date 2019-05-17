function simulateMonodomain(problem_name)
% This function simulates the monodomain equation on a mesh where diffusion
% coefficient is allowed to vary. No-flux regions may also be included
% inside the mesh. A vertex-centred approach is used for the actual
% integration step, as a method that deals more nicely with a non-diagonal
% diffusion coefficient.
%
% Input is the "problem_name", the filename of the problem file (.mat) that
% defines the problem to be simulated. These .mat files are created by
% separate functions (see example functions createOpenTissueProblem and
% createWedgeProblem to see how these are created)

% Monodomain parameters
lambda = 1;                               % Ratio of intracellular and extracellular conducitivies
chi = 2000;                               % Surface-to-volume ratio for tissue (cm^-1)
Cm = 1;                                   % Tissue capacitance per unit area (?F/cm²) - cell capacitance is defined in ionic model files

% Stimulus settings
stim_dur = 1;                             % Stimulus duration (ms)
stim_amp = 52;                            % Amplitude of stimulus per unit area (?A/cm²)
stim_times1 = [20];  % Vector of times to stimulate sites marked as a primary stimulus (ms)
stim_times2 = [];  % Vector of times to stimulate sites marked as a secondary stimulus (ms)

% Timestepping and solution methods
t_end = 1000;                             % Simulation time (ms)
dt = 0.01;                                % Timestep (ms)
reac_per_diffuse = 1;                     % Number of reaction steps to perform before each diffusive update
diff_exact = 0;                           % Require exact solves (direct methods) for linear system in diffusive updates

% Plotting
visualise = 1;                            % Flag for whether to visualise or not
save_anim = 1;                            % Flag for whether or not to save an animation (filename same as problem name, CAREFUL not to overwite!)
plot_interval = 2.5;                       % Time interval for plotting (ms)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the problem defined by the provided filename
load([problem_name,'.mat'], 'problem');

% Calculate the scale factor for diffusion (in monodomain model)
scale = lambda / (lambda + 1) / chi / Cm;

% Build the matrix representing the linear system constructed for diffusive
% updates - this code also outputs a list of which sites are actually
% active in the model
[A, b, active] = encodeDiffusiveProblem(problem.D_tensor, problem.Vfrac, problem.grid, dt, scale);
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

% Create video object if one is needed
if save_anim && visualise
    vid_obj = VideoWriter([problem_name,'.avi']);
    open(vid_obj);
end

% Initialise problem
[V, S] = initialiseProblem(cell_models, model_assignments, active);

% Initialise timers
diff_time = 0; reac_time = 0;

% Loop over time integrations
t = 0;
while t < t_end
    
    % Reaction step (perform a number of these for each diffusive update)
    for j = 1:reac_per_diffuse
    
        % Increment time
        t = t + dt;
        
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
        tic;
        [I_ion, S_active] = processReaction(V(active), S(active,:), dt, I_stim(active), cell_models, model_assignments(active));
        reac_time = reac_time + toc;
        
        % Update V and S at the active sites using the calculated value
        S(active,:) = S_active;
        V(active) = V(active) - dt * (1/Cm) * (I_ion + I_stim(active));
        
    end
        
    % Process diffusive update - uses current voltage values and the
    % diffusive update matrix A
    tic;
    V_active = processDiffusion(V(active), A, b, dt, diff_exact);
    V(active) = V_active;
    diff_time = diff_time + toc;
    
    % Check plot frequency, plot if hit
    if visualise && ( t - floor(t / plot_interval) * plot_interval <= dt*reac_per_diffuse )
       
        scatter(nodeX, nodeY, 25, V, 'filled' );
        caxis([-90 40]);
        title(['Voltage map at time t = ',num2str(t)], 'FontSize', 24);
        
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

% Output finish message, and running times
fprintf('Complete! Spent \n%g seconds on diffusion \n %g seconds on reaction \n \n', diff_time, reac_time);