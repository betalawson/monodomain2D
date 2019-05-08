function [V, S] = initialiseProblem(cell_models, model_assignments, active)
% Initialises the problem at a number of points specified by N, according
% to the model specified

% Loop over all cell models present in this simulation, processing all 
% cells assigned to a single model as one batch
V = zeros(size(model_assignments));
for k = 1:length(cell_models)
    
    % Find which cells use the current model
    batch = (model_assignments == k);
    % Perform a Rush-Larsen update using this model for these cells
    try
        [V(batch), S(batch,:)] = feval(['initialise',cell_models{k}], sum(batch));
    catch
        error('Initialisation failed for model %s.\nConfirm file %s.m exists, and verify it has no errors \n', cell_models{k}, ['initialise',cell_models{k}]);
    end
    
end


% Turn off all inactive sites (NaN is used for removal from plotting)
V(~active) = NaN;

