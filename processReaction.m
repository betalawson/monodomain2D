function [I_ion, S, Sinf, invtau] = processReaction(V, S, S_old, Sinf_old, invtau_old, dt, I_stim, I_stim_old, cell_models, model_assignments, mesh, second_order)
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
    
    % Only process if there are some nodes actually using this model
    if sum(batch) > 0
        
        try 
            % Perform the timestepping update for these nodes. Call a
            % different function depending on whether using second order
            % timestepping or not
            if second_order
                
                % If this is the first step (old information not present),
                % don't try to read out individual nodes from it. 
                % Otherwise, do.
                if isempty(Sinf_old)  
                    [I_ion(batch), S_batch, Sinf_batch, invtau_batch] = feval(['SecondOrderUpdate',this_model], V(batch), S(batch,:), S_old(batch,:), [], [], dt, I_stim(batch), I_stim_old(batch,:));
                else
                    [I_ion(batch), S_batch, Sinf_batch, invtau_batch] = feval(['SecondOrderUpdate',this_model], V(batch), S(batch,:), S_old(batch,:), Sinf_old(batch,:), invtau_old(batch,:), dt, I_stim(batch), I_stim_old(batch,:));
                end
                
            else
                [I_ion(batch), S_batch] = feval(['RLUpdate',this_model], V(batch), S(batch,:), dt, I_stim(batch));
            end
                
            % Store the batch outputs (performed this way so that cell
            % models with different numbers of state variables are
            % compatible
            S(batch,1:size(S_batch,2)) = S_batch;
            
            % Also update Sinf and invtau values for use in second order
            % timestepping
            if second_order
                Sinf(batch,1:size(Sinf_batch,2)) = Sinf_batch;
                invtau(batch,1:size(invtau_batch,2)) = invtau_batch;
            end
            
        catch 
            % Give a more meaningful error to user if the above failed
            if second_order
                func_prefix = 'SecondOrderUpdate';
            else
                func_prefix = 'RLUpdate';
            end
            error('Processing reaction term failed for model %s.\nConfirm file %s.m exists, and verify it has no errors \n', cell_models{k}, [func_prefix,cell_models{k}]);
        
        end
    end
end

end