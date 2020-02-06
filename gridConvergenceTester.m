function [mesh_separations, CVs, CVs_reactiononly, CVs_localonly] = gridConvergenceTester

% Set up the different gridsizes for which the propagation speed should be
% calculated
mesh_separations = [1000, 500];

% Initialise the vector of conduction velocities
CVs = zeros(length(mesh_separations),1);
CVs_reactiononly = zeros(length(mesh_separations),1);
CVs_localonly = zeros(length(mesh_separations),1);

% Loop over each different mesh separation
for k = 1:length(mesh_separations)
    
    % Create the problem, pause for safety, and then load it
    createOpenTissueProblem('test',mesh_separations(k));
    pause(2);
    load('test.mat','problem');
    
    % Run the wavespeed calculation code
    CVs(k) = calculateCV(problem);
    
    % Move into separate version of code that only nonlocally integrates
    % the reaction term, run that code, then return folder
    cd('ReactionOnly');
    CVs_reactiononly(k) = calculateCV(problem);
    cd('..');
    
    % Move into separate version of code that only nonlocally integrates
    % the reaction term, run that code, then return folder
    cd('LocalOnly');
    CVs_localonly(k) = calculateCV(problem);
    cd('..');
    
    fprintf('Completed run with mesh separation %d microns\n', mesh_separations(k));
    
end

