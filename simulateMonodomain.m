function simulateMonodomain
% This function simulates the monodomain equation on a mesh where diffusion
% coefficient is allowed to vary. No-flux regions may also be included
% inside the mesh. A vertex-centred approach is used for the actual
% integration step, as a method that deals more nicely with a non-diagonal
% diffusion coefficient

% Monodomain parameters
lambda = 1.4;                           % Ratio of intracellular and extracellular conducitivies
chi = 2000;                             % Surface-to-volume ratio for tissue (cm^-1)
Cm = 1;                                 % Tissue capacitance per unit area (?F/cm²) - cell capacitance is defined in ionic model files

% Problem specification
problem_name = 'leftstim';               % Filename for the problem to solve (create problem .mat files using the "create" functions)
model = 'TT04epi';                       % Ionic model

% Stimulus settings
stim_dur = 1;                           % Stimulus duration (ms)
stim_amp = 52;                          % Amplitude of stimulus per unit area (?A/cm²)
stim_times = [20, 520, 920, 1300, 1680];% Vector of times to stimulate (ms)

% Timestepping
t_end = 2000;                           % Simulation time (ms)
dt = 0.01;                              % Timestep (ms)
reac_per_diffuse = 1;                   % Number of reaction steps to perform before each diffusive update

% Plotting
visualise = 1;                          % Flag for whether to visualise or not
save_anim = 1;                          % Flag for whether or not to save an animation (filename same as problem name, CAREFUL not to overwite!)
plot_interval = 10;                     % Time interval for plotting (ms)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the problem defined by the provided filename
load([problem_name,'.mat'], 'problem');

% Build the matrix representing the linear system constructed for diffusive
% updates - this code also outputs a list of which sites are actually
% active in the model
[A, b, active] = encodeDiffusiveProblem(problem.D_tensor, problem.Vfrac, problem.grid);

% Calculate the scale factor for diffusion (in monodomain model)
scale = lambda / (lambda + 1) / chi / Cm;

% Read out the number of nodes to solve at
N = length(active);

% Initialise problem
[V, S] = initialiseProblem(N, model, active);

% Read out stimulus sites, and convert to a vectorised form
stim_sites = problem.stim_sites;
stim_sites = stim_sites';
stim_sites = stim_sites(:);

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
        if any( (t - stim_times) <= stim_dur & (t - stim_times) >= 0 )
            I_stim(stim_sites) = -stim_amp;
        end
    
        % Process reaction update - uses current voltage values and current
        % state variable values, S. Only processes active sites
        tic;
        [I_ion, S_active] = processReaction(V(active), S(active,:), dt, I_stim(active), model);
        reac_time = reac_time + toc;
        
        % Update V and S at the active sites using the calculated value
        S(active,:) = S_active;
        V(active) = V(active) - dt * (1/Cm) * (I_ion + I_stim(active));
        
    end
        
    % Process diffusive update - uses current voltage values and the
    % diffusive update matrix A
    tic;
    V_active = processDiffusion(V(active), A, b, dt, scale);
    V(active) = V_active;
    diff_time = diff_time + toc;
    
    % Check plot frequency, plot if hit
    if visualise && ( t - floor(t / plot_interval) * plot_interval <= dt*reac_per_diffuse )
       
        scatter(nodeX, nodeY, 10, V, 'filled' );
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