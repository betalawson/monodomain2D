function [I_ion, S] = processReaction(V, S, dt, I_stim, cell_models, model_assignments, mesh, parallel)
% This function applies Rush-Larsen integration to process a timestep of
% length dt, according to the cell model specified in input variable
% 'model'. The outputs are the new state variables, and the vector of total
% currents at each node,  I_ion.
%
% I_stim is provided simply so models that take it into account as a flow
% of potassium ions can accurately represent this.

% Initialise the current output vector
I_ion = zeros(size(model_assignments));

% Loop over all cell models present in this simulation, processing all
% cells assigned to a single model as one batch
for k = 1:length(cell_models)
    
    % Find which cells use the current model
    batch = (model_assignments == k & mesh.active);
    % Store current model as a variable (just to remove a warning)
    this_model = cell_models{k};
    
    % Perform a Rush-Larsen update using this model for these nodes
    % (assuming some nodes are marked with this model)
    if sum(batch) > 0
        
        if ~parallel % If not using parallel, process all at once
            try
                [I_ion(batch), S(batch,:)] = feval(['RLUpdate',this_model], V(batch), S(batch,:), dt, I_stim(batch));
            catch
                error('Processing reaction term failed for model %s.\nConfirm file %s.m exists, and verify it has no errors \n', cell_models{k}, ['RLUpdate',cell_models{k}]);
            end
            
        else  % Otherwise, split up the work onto multiple cores
            
            % Read out number of cores
            Ncores = feature('numcores');
            
            % Initialise job variables to remove a warning
            node_jobs = cell(Ncores,1);
            V_jobs = cell(Ncores,1);
            S_jobs = cell(Ncores,1);
            I_stim_jobs = cell(Ncores,1);
            
            % Divide workload into this many pieces
            process_cells = find(batch);
            Njobs = length(process_cells);
            for m = 1:Ncores
                node_jobs{m} = process_cells(m:Ncores:Njobs);
                V_jobs{m} = V(node_jobs{m});
                S_jobs{m} = S(node_jobs{m},:);
                I_stim_jobs{m} = I_stim(node_jobs{m});
            end
            
            % Perform jobs on each individual processor
            parfor m = 1:Ncores
                try
                    [I_results{m}, S_results{m}] = feval(['RLUpdate',this_model], V_jobs{m}, S_jobs{m}, dt, I_stim_jobs{m});
                catch
                    error('Processing reaction term failed for model %s.\nConfirm file %s.m exists, and verify it has no errors \n', this_model, ['RLUpdate',cell_models{k}]);
                end
            end
            
            % Recollect the results
            for m = 1:Ncores
                I_ion(node_jobs{m}) = I_results{m};
                S(node_jobs{m},:) = S_results{m};
            end
            
        end
    end
    
end

end