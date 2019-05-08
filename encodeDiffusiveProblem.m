function [A, b, active] = encodeDiffusiveProblem( D_tensor, Vfrac, grid, dt, scale )

% First, create a map informing which of the four elements surrounding
% every node is filled with fibrosis or not (causing it not to contribute
% to the integral of fluxes). No-flux boundaries around the edges are
% handled the exact same way, by creating a layer of 'fibrosis' around the
% whole problem domain. Because the diffusive term is:
%
%  ( 1 / Vfrac )
%
% with Vfrac the volume fraction, it is an input and factors into the
% equations calculated below. The value of (1 / Vfrac) at each node is
% calculated by averaging Vfrac and then taking the reciprocal of the
% average

% Specify here a flag that controls how volume fraction averages are
% calculated
include_occlusions = 0;     % 0 - do not include occlusions, only available space contributes to calculated values of volume fraction
                            % 1 - include internal occlusions, but disclude the no-flux boundaries from calculations
                            % 2 - include all occlusions, including no-flux boundary barriers, in calculations

% Read out individual diffusion tensor elements (matrices corresponding to
% their values throughout the domain)
D_xx = D_tensor.D_xx;
D_xy = D_tensor.D_xy;
D_yy = D_tensor.D_yy;

% Read out the grid sizes
dx = grid.dx;
dy = grid.dy;


% Read out the size of the input matrices
[Ny, Nx] = size(D_xx);

% Create an occupancy map that tracks which locations are actually
% non-diffusive (fibrotic)
occ_map = ( abs(D_xx) + abs(D_xy) + abs(D_yy) ) < 1e-10;

% Add a layer of fibrosis around the occupancy map to ease the
% implementation of no-flux boundary conditions on the outside
occ_map_ext = [ones(1,Nx+2); [ones(Ny,1), occ_map, ones(Ny,1)]; ones(1,Nx+2) ];

% Match this with a layer of zero diffusivity around the diffusivity maps
% provided
D_xx_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_xx, zeros(Ny,1)]; zeros(1,Nx+2) ];
D_xy_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_xy, zeros(Ny,1)]; zeros(1,Nx+2) ];
D_yy_ext = [zeros(1,Nx+2); [zeros(Ny,1), D_yy, zeros(Ny,1)]; zeros(1,Nx+2) ];

% Similarly, set the volume fractions as non-contributing in the case of
% the outer boundaries by marking boundaries as nan in a separate matrix,
% and giving them a volume fraction of zero
Vfrac_ext = [zeros(1,Nx+2); [zeros(Ny,1), Vfrac, zeros(Ny,1)]; zeros(1,Nx+2) ];
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
             
% Convert these matrices into a vector for ease of use in the matrix system
Vfrac_ext = Vfrac_ext'; Vfrac_ext = Vfrac_ext(:);
Vfrac_nodes = Vfrac_nodes'; Vfrac_nodes = Vfrac_nodes(:);
             
% Define the number of nodes in the normal and extended node grids
Nn = (Nx+1)*(Ny+1);
Nf = (Nx+3)*(Ny+3);

% Calculate the dependences of the bilinearly interpolated flux values at 
% the four control volume boundary midpoints for each non-fibrotic element.
% By dependencies, this refers to how these fluxes depend on each of the
% four surrounding note points. To do this, we use a vectorised approach
% (i.e. occupancy matrix becomes a vector)

% Convert occupancy and diffusivities to vectors
occ_v = occ_map_ext'; occ_v = occ_v(:);
D_xx_v = D_xx_ext'; D_xx_v = D_xx_v(:);
D_xy_v = D_xy_ext'; D_xy_v = D_xy_v(:);
D_yy_v = D_yy_ext'; D_yy_v = D_yy_v(:);

% Vertical fluxes through the control volume boundaries [dl, dr, ul, ur]
J_W = Vfrac_ext .* ~occ_v .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dy/2 - D_yy_v * 3*dx/4, D_xy_v * dy/2 - D_yy_v * dx/4, -D_xy_v * dy/2 + D_yy_v * 3*dx/4, D_xy_v * dy/2 + D_yy_v * dx/4 ] );
J_E = Vfrac_ext .* ~occ_v .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dy/2 - D_yy_v * dx/4, D_xy_v * dy/2 - D_yy_v * 3*dx/4, -D_xy_v * dy/2 + D_yy_v * dx/4, D_xy_v * dy/2 + D_yy_v * 3*dx/4 ] );
J_S = Vfrac_ext .* ~occ_v .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dx/2 - D_xx_v * 3*dy/4, -D_xy_v * dx/2 + D_xx_v * 3*dy/4, D_xy_v * dx/2 - D_xx_v * dy/4, D_xy_v * dx/2 + D_xx_v * dy/4 ] );
J_N = Vfrac_ext .* ~occ_v .* ( 1 / (dx*dy) ) .* ( [ -D_xy_v * dx/2 - D_xx_v * dy/4, -D_xy_v * dx/2 + D_xx_v * dy/4, D_xy_v * dx/2 - D_xx_v * 3*dy/4, D_xy_v * dx/2 + D_xx_v * 3*dy/4 ] );

% First create lists of the i and j co-ordinates of all *node* points
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

% Also create a remapped vector that converts the list of all node
% locations to their locations in a larger (Nx+3) by (Ny+3) grid of nodes
% (so that everything can be done with matrix indexing)
nlist = ( j * (Nx+3) + i + 1 );

% Now, the control volume formulation of the diffusive part of the problem
% is used to build the matrix, in terms of the above fluxes. The sparse
% matrix is created by first creating vectors [i,j,v] that define the row,
% column and value of the elements, respectively. This saves a great deal
% of memory and time
i_v = []; j_v = []; val_v = [];

% Node's dependence on itself
i_v = [i_v, nlist]; j_v = [j_v, nlist];
val_v = [val_v; ( dx/2 * J_W(eloc_ur,1) + dy/2 * J_S(eloc_ur,1) + dx/2 * J_E(eloc_ul,2) - dy/2 * J_S(eloc_ul,2) - dx/2 * J_W(eloc_dr,3) + dy/2 * J_N(eloc_dr,3) - dx/2 * J_E(eloc_dl,4) - dy/2 * J_N(eloc_dl,4) ) ./ Vfrac_nodes ];

% Node's dependence on North-East node
i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)+1];
val_v = [val_v; ( dx/2 * J_W(eloc_ur,4) + dy/2 * J_S(eloc_ur,4) ) ./ Vfrac_nodes];

% Node's dependence on North-West node
i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)-1];
val_v = [val_v; ( dx/2 * J_E(eloc_ul,3) - dy/2 * J_S(eloc_ul,3) ) ./ Vfrac_nodes];

% Node's dependence on South-East node
i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)+1];
val_v = [val_v; (- dx/2 * J_W(eloc_dr,2) + dy/2 * J_N(eloc_dr,2) ) ./ Vfrac_nodes];

% Node's dependence on South-West node
i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)-1];
val_v = [val_v; (- dx/2 * J_E(eloc_dl,1) - dy/2 * J_N(eloc_dl,1) ) ./ Vfrac_nodes];

% Node's dependence on East node
i_v = [i_v, nlist]; j_v = [j_v, nlist+1];
val_v = [val_v; (dx/2 * J_W(eloc_ur,2) + dy/2 * J_S(eloc_ur,2) - dx/2 * J_W(eloc_dr,4) + dy/2 * J_N(eloc_dr,4) ) ./ Vfrac_nodes];

% Node's dependence on West node
i_v = [i_v, nlist]; j_v = [j_v, nlist-1];
val_v = [val_v; (dx/2 * J_E(eloc_ul,1) - dy/2 * J_S(eloc_ul,1) - dx/2 * J_E(eloc_dl,3) - dy/2 * J_N(eloc_dl,3) ) ./ Vfrac_nodes];

% Node's dependence on North node
i_v = [i_v, nlist]; j_v = [j_v, nlist+(Nx+3)];
val_v = [val_v; (dx/2 * J_W(eloc_ur,3) + dy/2 * J_S(eloc_ur,3) + dx/2 * J_E(eloc_ul,4) - dy/2 * J_S(eloc_ul,4) ) ./ Vfrac_nodes];

% Node's dependence on South node
i_v = [i_v, nlist]; j_v = [j_v, nlist-(Nx+3)];
val_v = [val_v; (- dx/2 * J_W(eloc_dr,1) + dy/2 * J_N(eloc_dr,1) - dx/2 * J_E(eloc_dl,2) - dy/2 * J_N(eloc_dl,2) ) ./ Vfrac_nodes];

% Before creating the matrix, set any values that appear to be zero (but
% are not due to rounding error) to explicitly zero so they are not stored
%val_v( abs(val_v) < 100*eps ) = 0;

% Now create the matrix
A = sparse(i_v, j_v, val_v, Nf, Nf);

% Scale the whole matrix by the volume of the control volume
A = A / (dx * dy);

% Now, grab out only the rows of the matrix that correspond to real nodes
% (i.e. without the dummy boundaries)
A = A(nlist,nlist);

% No source term with straight diffusion, as implemented here
b = zeros(Nn,1);

% Track which nodes are 'active' (not buried completely in fibrosis)
active = ~all([occ_v(eloc_dl), occ_v(eloc_dr), occ_v(eloc_ul), occ_v(eloc_ur) ], 2);

% Reduce diffusive update matrix and vector according to these nodes. This
% also removes rows of zeros from the matrix, allowing use of backslash
A = A(active, active);
b = b(active);

% Finally, convert the matrix to the form actually used in solving the
% diffusive updates (performed once here to save time)
% That is, instead of the matrix defined above, which is for form
% dV/dt (diffusive) =  A V - b,
% it is now converted to form
% A V_new = V_old - dt * b.
% This is also the moment where the 'scale' of diffusion (factors in the
% monodomain equation) are introduced
A = speye(size(A)) - dt * scale * A;

end

