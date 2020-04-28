function createOnePatchProblem(filename)
% This function creates a single element with 4 nodes at its corners, that
% do not interact diffusively. This is used for single cell problems.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid separation is set automatically (irrelevant)
mesh_separation = 100;           % (in microns)

% Define the diffusion tensor
D = [ 0, 0; 0, 0 ];    % Zero diffusion for a 0D problem

% Convert the grid separation to centimetres for consistent units
mesh_separation = mesh_separation / 10000;

% Use one element
Nx = 1;
Ny = 1;
Lx = Nx * mesh_separation;
Ly = Ny * mesh_separation;

% Initialise occupancy map
occ_map = false;


%%% Define volumes of each element
Vfrac = double(~occ_map);      % Volume fraction of one in all non-occupied elements (no elements are partially occupied)


%%% Define node locations (element corners)
[nodeX, nodeY] = meshgrid( linspace(0,Lx,Nx+1), linspace(0,Ly,Ny+1) );

%%% All sites are stimulus sites
stim_sites1 = true(size(nodeY));
stim_sites2 = true(size(nodeY));


%%% Specify the cell model to use at all sites

% List cell models that will be used here
cell_models = {'CRN'};
% Assign models to cells (all nodes are model 1)
model_assignments = ones(size(nodeX));



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
problem.grid.dx = mesh_separation;
problem.grid.dy = mesh_separation;
problem.grid.Lx = Lx;
problem.grid.Ly = Ly;
problem.Nx = Nx;
problem.Ny = Ny;
nodeX = nodeX'; nodeX = nodeX(:);
nodeY = nodeY'; nodeY = nodeY(:);
problem.nodeX = nodeX;
problem.nodeY = nodeY;
problem.stim_sites1 = stim_sites1;
problem.stim_sites2 = stim_sites2;
problem.cell_models = cell_models;
problem.model_assignments = model_assignments;

% Save the problem
save([filename,'.mat'],'problem');

end

