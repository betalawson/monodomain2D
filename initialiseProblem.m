function [V, S] = initialiseProblem(N, model, active)
% Initialises the problem at a number of points specified by N, according
% to the model specified

% Initialise voltage and state variable values according to selected model
switch model
    
    case {'TT04epi', 'TT04endo', 'TT04M'}
        [V, S] = initialiseTT04(N);
        
    case {'TT06epi', 'TT06endo', 'TT06M'}
        [V, S] = initialiseTT06(N);
        
    case {'TT3epi', 'TT3endo', 'TT3M'}
        [V, S] = initialiseTT3(N);

end

% Turn off all inactive sites
V(~active) = NaN;

