function [K, M, mesh] = encodeProblem(occ_map, D_tensor, Vfrac, grid, alpha)
% This function calculates the stiffness and mass matrices associated with
% the input information regarding the mesh (including what sites are
% occupied, and the diffusivities and volume fractions)

% Specify whether to approximate a nodal volume fraction or to incorporate
% the element-associated volume fractions into the mass matrix
volume_weighted_mass_matrix = 1;

% Specify here a flag that controls how volume fraction averages are
% calculated. Only used when nodal averages are required (if the above
% setting is turned off)
%     0 - do not include occlusions, only available space contributes to calculated values of volume fraction
%     1 - include internal occlusions, but disclude the no-flux boundaries from calculations
%     2 - include all occlusions, including no-flux boundary barriers, in calculations
include_occlusions = 0;



%%% LOTS OF SETUP

% Read out the individual diffusion tensor elements (matrices corresponding
% to their values throughout the domain)
D_xx = D_tensor.D_xx;
D_xy = D_tensor.D_xy;
D_yy = D_tensor.D_yy;

% Read out the grid sizes
dx = grid.dx; dy = grid.dy;

% Read out the size of the input matrices
[Ny, Nx] = size(D_xx);

% In order to ease implementation of no-flux boundaries, add a layer of
% "occlusion" around the occupancy map
occ_map_ext = [true(1,Nx+2); [true(Ny,1), occ_map, true(Ny,1)]; true(1,Nx+2) ];

% Match this with a layer of zero diffusivity around the diffusivity maps
% provided
D_xx_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_xx, zeros(Ny,1)]; zeros(1,Nx+2) ];
D_xy_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_xy, zeros(Ny,1)]; zeros(1,Nx+2) ];
D_yy_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_yy, zeros(Ny,1)]; zeros(1,Nx+2) ];

% Create an extended matrix of the volume fraction (a layer of zeros around
% edges to represent outside the domain as fully occluded)
Vfrac_ext = [zeros(1,Nx+2); [zeros(Ny,1), Vfrac, zeros(Ny,1)]; zeros(1,Nx+2) ];

% A very similar calculation calculates the volume of each control volume
CV_vols = dx * dy * ( ~occ_map_ext(1:Ny+1, 1:Nx+1) + ~occ_map_ext(2:Ny+2, 1:Nx+1) + ~occ_map_ext(1:Ny+1, 2:Nx+2) + ~occ_map_ext(2:Ny+2, 2:Nx+2) ) / 4;

% Convert these matrices into a vector for ease of use in the matrix system
Vfrac_ext = Vfrac_ext'; Vfrac_ext = Vfrac_ext(:);
CV_vols = CV_vols'; CV_vols = CV_vols(:);
% Also make vector versions of the occupancy and diffusivities
occ_v_ext = occ_map_ext'; occ_v_ext = occ_v_ext(:);
D_xx_v = D_xx_ext'; D_xx_v = D_xx_v(:);
D_xy_v = D_xy_ext'; D_xy_v = D_xy_v(:);
D_yy_v = D_yy_ext'; D_yy_v = D_yy_v(:);

% Define the number of nodes in the normal and extended node grids
Nn = (Nx+1)*(Ny+1);
Nf = (Nx+3)*(Ny+3);

% Create lists of the i and j co-ordinates of all *node* points
j = ceil( (1:Nn) / (Nx+1) );
i = ( (1:Nn) - (j-1)*(Nx+1) );

% Create a vector of the corresponding position in the *element* grid
% (which is Ny+2 by Nx+2 after extension with boundaries). This position
% refers to the element down and left from the node point
eloc_dl = ( (j-1) * (Nx+2) + i );

% Calculate the corresponding positions of other elements accordingly
eloc_dr = eloc_dl + 1;
eloc_ul = eloc_dl + (Nx+2);
eloc_ur = eloc_dl + (Nx+2) + 1;

% Track which nodes are 'active' (not buried completely in fibrosis), and
% 'free' (no fully-occluded elements present)
active = ~all([occ_v_ext(eloc_dl), occ_v_ext(eloc_dr), occ_v_ext(eloc_ul), occ_v_ext(eloc_ur) ], 2);

% Get lists of free and non-free nodes
active_nodes = find(active);



%%% MASS MATRIX

% Use the occupancy of elements to construct the mass matrix. The mass
% matrix here uses a piecewise bilinear interpolation
i_v = []; j_v = []; val_v = [];


% Each piecewise bilinear interpolation contributes a mass matrix stencil 
% of form:
%       [   0      3/64    1/64  ]
%       [   0      9/64    3/64  ]
%       [   0       0       0    ]
% If volume fractions are included in mass matrix calculations, each piece
% of this form will be weighted by the volume fraction, otherwise they are
% weighted by ~occ (which simply denotes whether there is a contribution or
% not)
if volume_weighted_mass_matrix
    
    % Centre node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes];
    val_v = [val_v; 9/64 * (Vfrac_ext(eloc_dl(active_nodes)) + Vfrac_ext(eloc_ul(active_nodes)) + Vfrac_ext(eloc_dr(active_nodes)) + Vfrac_ext(eloc_ur(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Right node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+1];
    val_v = [val_v; 3/64 * (Vfrac_ext(eloc_dr(active_nodes)) + Vfrac_ext(eloc_ur(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Left node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-1];
    val_v = [val_v; 3/64 * (Vfrac_ext(eloc_dl(active_nodes)) + Vfrac_ext(eloc_ul(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Up node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+(Nx+1)];
    val_v = [val_v; 3/64 * (Vfrac_ext(eloc_ul(active_nodes)) + Vfrac_ext(eloc_ur(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Down node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-(Nx+1)];
    val_v = [val_v; 3/64 * (Vfrac_ext(eloc_dl(active_nodes)) + Vfrac_ext(eloc_dr(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Up-right node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+1+(Nx+1)];
    val_v = [val_v; 1/64 * Vfrac_ext(eloc_ur(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Up-left node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-1+(Nx+1)];
    val_v = [val_v; 1/64 * Vfrac_ext(eloc_ul(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Down-right node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+1-(Nx+1)];
    val_v = [val_v; 1/64 * Vfrac_ext(eloc_dr(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Down-left node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-1-(Nx+1)];
    val_v = [val_v; 1/64 * Vfrac_ext(eloc_dl(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
else
    
    % Centre node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes];
    val_v = [val_v; 9/64 * (~occ_v_ext(eloc_dl(active_nodes)) + ~occ_v_ext(eloc_ul(active_nodes)) + ~occ_v_ext(eloc_dr(active_nodes)) + ~occ_v_ext(eloc_ur(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Right node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+1];
    val_v = [val_v; 3/64 * (~occ_v_ext(eloc_dr(active_nodes)) + ~occ_v_ext(eloc_ur(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Left node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-1];
    val_v = [val_v; 3/64 * (~occ_v_ext(eloc_dl(active_nodes)) + ~occ_v_ext(eloc_ul(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Up node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+(Nx+1)];
    val_v = [val_v; 3/64 * (~occ_v_ext(eloc_ul(active_nodes)) + ~occ_v_ext(eloc_ur(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Down node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-(Nx+1)];
    val_v = [val_v; 3/64 * (~occ_v_ext(eloc_dl(active_nodes)) + ~occ_v_ext(eloc_dr(active_nodes)) ) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Up-right node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+1+(Nx+1)];
    val_v = [val_v; 1/64 * ~occ_v_ext(eloc_ur(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Up-left node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-1+(Nx+1)];
    val_v = [val_v; 1/64 * ~occ_v_ext(eloc_ul(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Down-right node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes+1-(Nx+1)];
    val_v = [val_v; 1/64 * ~occ_v_ext(eloc_dr(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
    % Down-left node
    i_v = [i_v; active_nodes];
    j_v = [j_v; active_nodes-1-(Nx+1)];
    val_v = [val_v; 1/64 * ~occ_v_ext(eloc_dl(active_nodes)) * dx * dy ./ CV_vols(active_nodes) ];
    
end

% Before constructing the sparse matrix, remove all terms from the vectors
% that have a zero value (some of these will reference locations outside
% the node grid, but it's fine to remove them because they are zero anyway)
retain_v = (val_v ~= 0);
i_v = i_v(retain_v);
j_v = j_v(retain_v);
val_v = val_v(retain_v);

% With the mass matrix fully formed, create it as a sparse matrix
M = sparse(i_v, j_v, val_v, Nn, Nn);

% Only grab out the active nodes
M = M(active, active);

%%% STIFFNESS MATRIX

% Calculate the dependences of the bilinearly interpolated flux values at
% the four control volume boundary midpoints for each non-fibrotic element.
% By dependencies, this refers to how these fluxes depend on each of the
% four surrounding note points. To do this, we use a vectorised approach
% (i.e. occupancy matrix becomes a vector)

% First create a remapped vector that converts the list of all node
% locations to their locations in a larger (Nx+3) by (Ny+3) grid of nodes
% (so that everything can be done with matrix indexing)
nlist = ( j * (Nx+3) + i + 1 );

% Vertical fluxes through the control volume boundaries [dl, dr, ul, ur]
J_W = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dy/2 - D_yy_v * 3*dx/4, D_xy_v * dy/2 - D_yy_v * dx/4, -D_xy_v * dy/2 + D_yy_v * 3*dx/4, D_xy_v * dy/2 + D_yy_v * dx/4 ] );
J_E = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dy/2 - D_yy_v * dx/4, D_xy_v * dy/2 - D_yy_v * 3*dx/4, -D_xy_v * dy/2 + D_yy_v * dx/4, D_xy_v * dy/2 + D_yy_v * 3*dx/4 ] );
J_S = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dx/2 - D_xx_v * 3*dy/4, -D_xy_v * dx/2 + D_xx_v * 3*dy/4, D_xy_v * dx/2 - D_xx_v * dy/4, D_xy_v * dx/2 + D_xx_v * dy/4 ] );
J_N = Vfrac_ext .* ~occ_v_ext .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dx/2 - D_xx_v * dy/4, -D_xy_v * dx/2 + D_xx_v * dy/4, D_xy_v * dx/2 - D_xx_v * 3*dy/4, D_xy_v * dx/2 + D_xx_v * 3*dy/4 ] );

% Initialise a new sparse matrix structure
i_v = []; j_v = []; val_v = [];



% Nodal volume fractions will or will not be included in the stiffness
% matrix depending on the specified setting
if volume_weighted_mass_matrix
    
    % Node's dependence on itself
    i_v = [i_v, nlist]; j_v = [j_v, nlist];
    val_v = [val_v; ( dx/2 * J_W(eloc_ur,1) + dy/2 * J_S(eloc_ur,1) + dx/2 * J_E(eloc_ul,2) - dy/2 * J_S(eloc_ul,2) - dx/2 * J_W(eloc_dr,3) + dy/2 * J_N(eloc_dr,3) - dx/2 * J_E(eloc_dl,4) - dy/2 * J_N(eloc_dl,4) ) ./ CV_vols ];
    
    % Node's dependence on North-East node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)+1];
    val_v = [val_v; ( dx/2 * J_W(eloc_ur,4) + dy/2 * J_S(eloc_ur,4) ) ./ CV_vols];
    
    % Node's dependence on North-West node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)-1];
    val_v = [val_v; ( dx/2 * J_E(eloc_ul,3) - dy/2 * J_S(eloc_ul,3) ) ./ CV_vols];
    
    % Node's dependence on South-East node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)+1];
    val_v = [val_v; (- dx/2 * J_W(eloc_dr,2) + dy/2 * J_N(eloc_dr,2) ) ./ CV_vols];
    
    % Node's dependence on South-West node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)-1];
    val_v = [val_v; (- dx/2 * J_E(eloc_dl,1) - dy/2 * J_N(eloc_dl,1) ) ./ CV_vols];
    
    % Node's dependence on East node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+1];
    val_v = [val_v; (dx/2 * J_W(eloc_ur,2) + dy/2 * J_S(eloc_ur,2) - dx/2 * J_W(eloc_dr,4) + dy/2 * J_N(eloc_dr,4) ) ./ CV_vols];
    
    % Node's dependence on West node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-1];
    val_v = [val_v; (dx/2 * J_E(eloc_ul,1) - dy/2 * J_S(eloc_ul,1) - dx/2 * J_E(eloc_dl,3) - dy/2 * J_N(eloc_dl,3) ) ./ CV_vols];
    
    % Node's dependence on North node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)];
    val_v = [val_v; (dx/2 * J_W(eloc_ur,3) + dy/2 * J_S(eloc_ur,3) + dx/2 * J_E(eloc_ul,4) - dy/2 * J_S(eloc_ul,4) ) ./ CV_vols];
    
    % Node's dependence on South node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)];
    val_v = [val_v; (- dx/2 * J_W(eloc_dr,1) + dy/2 * J_N(eloc_dr,1) - dx/2 * J_E(eloc_dl,2) - dy/2 * J_N(eloc_dl,2) ) ./ CV_vols];
    
else
    
    % Set the volume fractions as non-contributing in the case of
    % the outer boundaries by marking boundaries as nan in a separate matrix,
    % and giving them a volume fraction of zero
    switch include_occlusions
        case 0   % Disclude all occlusions
            include = logical(Vfrac_ext);
        case 1   % Only disclude boundaries
            include = [zeros(1,Nx+2); [zeros(Ny,1), ones(size(Vfrac)), zeros(Ny,1)]; zeros(1,Nx+2) ];
        case 2   % Include everything
            include = ones(size(Vfrac_ext));
    end
    
    % Now use these values to calculate the average volume fraction associated
    % with each node
    Vfrac_nodes = ( ...
        include(1:Ny+1, 1:Nx+1) .* Vfrac_ext(1:Ny+1, 1:Nx+1) ...
        + include(2:Ny+2, 1:Nx+1) .* Vfrac_ext(2:Ny+2, 1:Nx+1) ...
        + include(1:Ny+1, 2:Nx+2) .* Vfrac_ext(1:Ny+1, 2:Nx+2) ...
        + include(2:Ny+2, 2:Nx+2) .* Vfrac_ext(2:Ny+2, 2:Nx+2) ...
        ) ./ ( include(1:Ny+1, 1:Nx+1) + include(2:Ny+2, 1:Nx+1) + include(1:Ny+1, 2:Nx+2) + include(2:Ny+2, 2:Nx+2) );
    
    % Convert this into a vector for use in below calculations
    Vfrac_nodes = Vfrac_nodes'; Vfrac_nodes = Vfrac_nodes(:);
    
    % Node's dependence on itself
    i_v = [i_v, nlist]; j_v = [j_v, nlist];
    val_v = [val_v; ( dx/2 * J_W(eloc_ur,1) + dy/2 * J_S(eloc_ur,1) + dx/2 * J_E(eloc_ul,2) - dy/2 * J_S(eloc_ul,2) - dx/2 * J_W(eloc_dr,3) + dy/2 * J_N(eloc_dr,3) - dx/2 * J_E(eloc_dl,4) - dy/2 * J_N(eloc_dl,4) ) ./ (CV_vols .* Vfrac_nodes) ];
    
    % Node's dependence on North-East node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)+1];
    val_v = [val_v; ( dx/2 * J_W(eloc_ur,4) + dy/2 * J_S(eloc_ur,4) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on North-West node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)-1];
    val_v = [val_v; ( dx/2 * J_E(eloc_ul,3) - dy/2 * J_S(eloc_ul,3) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on South-East node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)+1];
    val_v = [val_v; (- dx/2 * J_W(eloc_dr,2) + dy/2 * J_N(eloc_dr,2) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on South-West node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)-1];
    val_v = [val_v; (- dx/2 * J_E(eloc_dl,1) - dy/2 * J_N(eloc_dl,1) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on East node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+1];
    val_v = [val_v; (dx/2 * J_W(eloc_ur,2) + dy/2 * J_S(eloc_ur,2) - dx/2 * J_W(eloc_dr,4) + dy/2 * J_N(eloc_dr,4) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on West node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-1];
    val_v = [val_v; (dx/2 * J_E(eloc_ul,1) - dy/2 * J_S(eloc_ul,1) - dx/2 * J_E(eloc_dl,3) - dy/2 * J_N(eloc_dl,3) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on North node
    i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)];
    val_v = [val_v; (dx/2 * J_W(eloc_ur,3) + dy/2 * J_S(eloc_ur,3) + dx/2 * J_E(eloc_ul,4) - dy/2 * J_S(eloc_ul,4) ) ./ (CV_vols .* Vfrac_nodes)];
    
    % Node's dependence on South node
    i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)];
    val_v = [val_v; (- dx/2 * J_W(eloc_dr,1) + dy/2 * J_N(eloc_dr,1) - dx/2 * J_E(eloc_dl,2) - dy/2 * J_N(eloc_dl,2) ) ./ (CV_vols .* Vfrac_nodes)];
    
end

% Before creating the matrix, set any values that appear to be zero (but
% are not due to rounding error) to explicitly zero so they are not stored
val_v( abs(val_v) < 100*eps ) = 0;

% Now create the matrix (and multiply by the scaling coefficient alpha)
K = alpha * sparse(i_v, j_v, val_v, Nf, Nf);

% Now, grab out only the rows of the matrix that correspond to real nodes
% (i.e. without the dummy boundaries)
K = K(nlist,nlist);

% And actually, we only want the active nodes
K = K(active, active);

%%% STORE MESH INFORMATION

mesh.dx = dx;
mesh.dy = dy;
mesh.Nx = Nx;
mesh.Ny = Ny;
mesh.occ_v_ext = occ_v_ext;
mesh.active = active;
mesh.CV_vols = CV_vols;

end