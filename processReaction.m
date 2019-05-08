function [I_ion, S] = processReaction(V, S, dt, I_stim, cell_models, model_assignments)
% This function applies Rush-Larsen integration to process a timestep of
% length dt, according to the cell model specified in input variable
% 'model'. 

% Loop over all cell models present in this simulation, processing all 
% cells assigned to a single model as one batch
I_ion = zeros(size(model_assignments));
for k = 1:length(cell_models)
    
    % Find which cells use the current model
    batch = (model_assignments == k);
    % Perform a Rush-Larsen update using this model for these cells
    try
        [I_ion(batch), S(batch,:)] = feval(['RLUpdate',cell_models{k}], V(batch), S(batch,:), dt, I_stim(batch));
    catch
        error('Processing reaction term failed for model %s.\nConfirm file %s.m exists, and verify it has no errors \n', cell_models{k}, ['RLUpdate',cell_models{k}]);
    end
    
end

end

