function [I_ion, S] = processReaction(V, S, dt, I_stim, cell_models, model_assignments, local, mesh)
% This function applies Rush-Larsen integration to process a timestep of
% length dt, according to the cell model specified in input variable
% 'model'. There is the option for nonlocal integration (see documentation
% and/or comments below), thus requiring the input of occupancy information

% Initialise current output, also the 'node' current output (these are
% different in the nonlocal case)
I_ion = zeros(size(model_assignments));
I_ion_nodes = zeros(size(model_assignments));

% Loop over all cell models present in this simulation, processing all 
% cells assigned to a single model as one batch
for k = 1:length(cell_models)
    
    % Find which cells use the current model
    batch = (model_assignments == k & mesh.active);
    % Perform a Rush-Larsen update using this model for these nodes
    % (assuming some nodes are marked with this model)
    if sum(batch) > 0
        try
            [I_ion_nodes(batch), S(batch,:)] = feval(['RLUpdate',cell_models{k}], V(batch), S(batch,:), dt, I_stim(batch));
        catch
            error('Processing reaction term failed for model %s.\nConfirm file %s.m exists, and verify it has no errors \n', cell_models{k}, ['RLUpdate',cell_models{k}]);
        end
    end
    
end


% Now, according to how to handle the integration of the reaction term over
% the control volume, calculate I_ion values that actually update the nodes
if local
    
    % Local integration just means using the value of I_ion at the CV
    % centre as the average value
    I_ion = I_ion_nodes;

else
    
    % Non-local integration means creating the bilinear interpolant between
    % I_ion values at the nodes, and integrating this over the CV. Only
    % non-occupied elements contribute to the integration
    
    % Read out mesh information
    dx = mesh.dx;
    dy = mesh.dy;
    Nx = mesh.Nx;
    Ny = mesh.Ny;
    active = mesh.active;
    
    % Find list of elements that must be processed (non-occupied sites)
    process_eles = find(~mesh.occ_v_ext);
    
    % Create a list of the nodes corresponding to each of these elements
    dl = process_eles - (Nx+2) - floor( process_eles / (Nx + 2) );
    dr = dl + 1;
    ul = dl + Nx+1;
    ur = dl + Nx+1 + 1;
    
    % Calculate the bilinear interpolant for each of these, expressed as
    %  F(:,1) + F(:,2) x + F(:,3) y + F(:,4) x y
    F = zeros( (Nx+2) * (Ny+2), 4);
    F(process_eles,:) = [ I_ion_nodes(dl),  ( I_ion_nodes(dr) - I_ion_nodes(dl) ) / dx, ( I_ion_nodes(ul) - I_ion_nodes(dl) ) / dy, ( I_ion_nodes(ur) + I_ion_nodes(dl) - I_ion_nodes(dr) - I_ion_nodes(ul) ) / (dx * dy) ];
    
    % Create a list of the elements corresponding to each active node
    active_nodes = find(active);
    dle = active_nodes + floor( (active_nodes - 1) / (Nx + 1) );
    dre = dle + 1;
    ule = dle + Nx+2;
    ure = dle + Nx+2 + 1;
    
    % Calculate I_ion as the integral (divided by CV volume) of these
    % interpolants over the different elements
    I_ion(active) = ( BilinInt(F(dle,:),dx,dy,1,1) + BilinInt(F(dre,:),dx,dy,0,1) + BilinInt(F(ule,:),dx,dy,1,0) + BilinInt(F(ure,:),dx,dy,0,0) ) ./ mesh.CV_vols;
    
end





function val = BilinInt( F, dx, dy, xplus, yplus )
% This function integrates the bilinear interpolations specified in F (each
% row of F a new interpolant) over a single quadrant (corresponding to a
% portion of the CV in a vertex-centred FVM setup). xplus and yplus specify
% the quadrant, values of 1 for both corresponding to top-right, and so on
val = F(:,1) * dx * dy / 4 + F(:,2) * (1 + 2*xplus) * dx^2 * dy / 16 + F(:,3) * (1 + 2*yplus) * dx * dy^2 / 16 + F(:,4) * (1 + 2*xplus) * (1 + 2*yplus) * dx^2 * dy^2 / 64;