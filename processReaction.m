function [I_ion, S, Sinf, invtau, b] = processReaction(V, S, Sinf_old, invtau_old, b_old, dt, I_stim, I_stim_old, cell_models, model_assignments, mesh, second_order, extra_params)
% This function applies Rush-Larsen integration to process a timestep of
% length dt, according to the cell model specified in input variable
% 'model'. The outputs are the new state variables, and the vector of total
% currents at each node,  I_ion.
%
% I_stim is provided simply so models that take it into account as a flow
% of potassium ions can accurately represent this.

% Initialise all output data
I_ion = zeros(size(model_assignments));

% Initialise output 

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
                    [I_ion(batch), S_batch, Sinf_batch, invtau_batch, b_batch] = feval(['SecondOrderUpdate',this_model], V(batch), S(batch,:), [], [], [], dt, I_stim(batch), I_stim_old(batch,:), extra_params);
                else
                    [I_ion(batch), S_batch, Sinf_batch, invtau_batch, b_batch] = feval(['SecondOrderUpdate',this_model], V(batch), S(batch,:), Sinf_old(batch,:), invtau_old(batch,:), b_old(batch,:), dt, I_stim(batch), I_stim_old(batch,:), extra_params);
                end
                
            else
                [I_ion(batch), S_batch] = feval(['RLUpdate',this_model], V(batch), S(batch,:), dt, I_stim(batch), extra_params);
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
                b(batch,1:size(b_batch,2)) = b_batch;
            else   % Otherwise set them to dummy stuff
                Sinf = NaN;
                invtau = NaN;
                b = NaN;
            end
            
        catch
            
            % Give a more meaningful error to user if the above failed
            if second_order
                func_prefix = 'SecondOrderUpdate';
            else
                func_prefix = 'RLUpdate';
            end
            fprintf('\n -------------- \n  FATAL ERROR: Processing reaction term failed for model %s.\n -------------- \n\nConfirm file %s.m exists, and verify it has no errors. \n \n Retrying call to indicate nature of error... \n\n', cell_models{k}, [func_prefix,cell_models{k}]);
        
            % Re-run model as called above so it outputs the error
            if second_order
                if isempty(Sinf_old) 
                    [I_ion(batch), S_batch, Sinf_batch, invtau_batch, b_batch] = feval(['SecondOrderUpdate',this_model], V(batch), S(batch,:), [], [], [], dt, I_stim(batch), I_stim_old(batch,:), extra_params);
                else
                    [I_ion(batch), S_batch, Sinf_batch, invtau_batch, b_batch] = feval(['SecondOrderUpdate',this_model], V(batch), S(batch,:), Sinf_old(batch,:), invtau_old(batch,:), b_old(batch,:), dt, I_stim(batch), I_stim_old(batch,:), extra_params);
                end
            else
                [I_ion(batch), S_batch] = feval(['RLUpdate',this_model], V(batch), S(batch,:), dt, I_stim(batch), extra_params);
            end
        
       end
    end
end

end