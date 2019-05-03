function [I_ion, S] = processReaction(V, S, dt, I_stim, model)
% This function applies Rush-Larsen integration to process a timestep of
% length dt, according to the cell model specified in input variable
% 'model'. 

switch model
    
    case 'TT04epi'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT04(V, S, dt, I_stim, 'epi');
    case 'TT04M'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT04(V, S, dt, I_stim, 'M');
    case 'TT04endo'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT04(V, S, dt, I_stim, 'endo');
    case 'TT06epi'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT06(V, S, dt, I_stim, 'epi');
    case 'TT06M'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT06(V, S, dt, I_stim, 'M');
    case 'TT06endo'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT06(V, S, dt, I_stim, 'endo');
    case 'TT3epi'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT3(V, S, dt, 'epi');
    case 'TT3M'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT3(V, S, dt, 'M');
    case 'TT3endo'    % Reduced Ten-Tusscher model (epicardial)
        [I_ion, S] = RLUpdateTT3(V, S, dt, 'endo');
        
    otherwise
        error('Model not found!');
        
end

end

