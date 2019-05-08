function createOpenTissueProblem(filename)
% This function creates a wedge problem

% Define the diffusion tensor
D = [ 1, 0; 0, 3 ];

% Define the number of elements in the problem
Nx = 50;
Ny = 50;

% Define the physical size of the problem (in centimetres)
Lx = 1;
Ly = 1;

% Specify stimulus region size (bottom left corner)
stim_radius = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Create a domain unoccupied by fibrosis. Volume fraction is one
%%% everywhere, occupancy (of blockage) zero everywhere

% Determine x and y co-ordinates of all element centres
dx = Lx / Nx;
dy = Ly / Ny;
[X,Y] = meshgrid( linspace(0,Lx-dx,Nx) + dx/2,   linspace(0,Ly-dy,Ny) + dy/2 );

% Initialise occupancy map
occ_map = zeros(size(Y));

%%% Define volumes of each element
Vfrac = ~occ_map;      % Volume fraction of one in all non-occupied elements (no elements are partially occupied)


%%% Create stimulus sites

% Nodes will be placed at element boundaries (vertex-centred finite volume)
% So first create node positions
[nodeX, nodeY] = meshgrid( linspace(0,Lx,Nx+1), linspace(0,Ly,Ny+1) );

% Initialise stimulus matrix to zeroes
stim_sites = false(size(nodeY));

% Set the bottom-left corner to be stimulated
stim_sites( sqrt(nodeX.^2 + nodeY.^2) <= stim_radius^2) = true;



%%% Specify the cell model to use at all sites

% List cell models that will be used here
cell_models = {'TT3epi'};
% Assign models to cells (by number)
model_assignments = zeros(size(nodeX));
model_assignments(:) = 1;




%%% Process and save all data

% Read out base diffusivity levels from the diffusion tensor
D_xx = D(1,1);
D_xy = D(1,2);
D_yy = D(2,2);

% Create matrices of diffusion values, with zero in blocked regions
D_xx = D_xx * (~occ_map);
D_xy = D_xy * (~occ_map);
D_yy = D_yy * (~occ_map);

% Store problem details in the 'problem' structure
problem.D_tensor.D_xx = D_xx;
problem.D_tensor.D_xy = D_xy;
problem.D_tensor.D_yy = D_yy;
problem.Vfrac = Vfrac;
problem.grid.dx = dx;
problem.grid.dy = dy;
problem.stim_sites = stim_sites;
problem.nodeX = nodeX;
problem.nodeY = nodeY;
problem.cell_models = cell_models;
problem.model_assignments = model_assignments;

% Save the problem
save([filename,'.mat'],'problem');




end

