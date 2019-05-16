function createWedgeProblem(filename)
% This function creates a wedge problem

% Define the diffusion tensor
D = [ 1, 0; 0, 3 ];

% Define the number of elements in the problem
Nx = 200;
Ny = 200;

% Define the physical size of the problem (in centimetres)
Lx = 1;
Ly = 1;

% Specify the horizontal location where the channel begins and ends
channel_start = 0;
channel_end = 0.5;
channel_constant_length = 0.5;
channel_start_width = 0.02;
channel_end_width = 0.02;

% Stimulus regions are the left and right edges of the domain (left
% primary, right secondary)
stim_width = 0.2;            % Width of stimulus regions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Create the pattern of the wedge-shaped blockage

% Determine x and y co-ordinates of all element centres
dx = Lx / Nx;
dy = Ly / Ny;
[X,Y] = meshgrid( linspace(0,Lx-dx,Nx) + dx/2,   linspace(0,Ly-dy,Ny) + dy/2 );

% Initialise occupancy map
occ_map = zeros(size(Y));

% Channel entrance
occ_map( X > channel_start & X < channel_start+channel_constant_length &  abs( Y - Ly/2 ) >= channel_start_width/2 ) = 1;

% Channel exit
occ_map( X < channel_end & X > channel_end-channel_constant_length &  abs( Y - Ly/2 ) >= channel_end_width/2 ) = 1;

% Channel opening region
occ_map( X >= channel_start+channel_constant_length & X <= channel_end-channel_constant_length & abs( Y - Ly/2) >= ( channel_start_width/2 + (channel_end_width - channel_start_width) / ( 2 * ( channel_end - channel_start - 2*channel_constant_length) ) * (X - channel_start - channel_constant_length) ) ) = 1;


%%% Define volumes of each element
Vfrac = ~occ_map;      % Volume fraction of one in all non-occupied elements (no elements are partially occupied)


%%% Create stimulus sites

% Nodes will be placed at element boundaries (vertex-centred finite volume)
% So first create node positions
[nodeX, nodeY] = meshgrid( linspace(0,Lx,Nx+1), linspace(0,Ly,Ny+1) );

% Initialise stimulus matrix to zeroes
stim_sites1 = false(size(nodeY));
stim_sites2 = false(size(nodeY));

% Set edges to be stimulus sites as requested
stim_sites1(nodeX <= stim_width) = true;
stim_sites2(nodeX >= Lx - stim_width) = true;



%%% Specify the cell model to use at all sites

% List cell models that will be used here
cell_models = {'TT3weak'};
% Assign models to cells (by number)
model_assignments = zeros(size(nodeX));
model_assignments(nodeX < 1) = 1;
model_assignments(nodeX >= 1) = 1;




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
problem.occ_map = occ_map;
problem.D_tensor.D_xx = D_xx;
problem.D_tensor.D_xy = D_xy;
problem.D_tensor.D_yy = D_yy;
problem.Vfrac = Vfrac;
problem.grid.dx = dx;
problem.grid.dy = dy;
problem.grid.Lx = Lx;
problem.grid.Ly = Ly;
problem.Nx = Nx;
problem.Ny = Ny;
problem.nodeX = nodeX;
problem.nodeY = nodeY;
problem.stim_sites1 = stim_sites1;
problem.stim_sites2 = stim_sites2;
problem.cell_models = cell_models;
problem.cell_models = cell_models;
problem.model_assignments = model_assignments;

% Save the problem
save([filename,'.mat'],'problem');




end

