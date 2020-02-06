# monodomain2D

### Code Features
A MATLAB implementation that solves the Monodomain equation in 2D. The code assumes a regular grid, but allows for the user to specify
heterogeneous diffusion tensors, and blocked regions. Different cell models can be specified for different nodes, but only the reduced
ten Tusscher *et al.* model (TT3 model) is included in the current implementation.
The monodomain equation solved is one that allows for "partial occlusion" by considering the intrinsic average of the typical monodomain 
equation, leading to a volume fraction term in the equation. For a standard monodomain implementation, simply use a volume fraction of unity
everywhere (and zero for occluded sites). The example problem files do this.

### Numerical Approach
The code uses a control volume finite element method (CV-FEM), in which a vertex-centred finite volume approach is used, but with interpolants
used to better approximate integrals over the control volume, giving an FEM-like quality. 

### Usage
Code is commented, and the majority of information regarding use of the code is provided via these comments. Basic procedure for usage is to
use the "create" file(s) as a template to see how to create your own problems in the structure required by the simulator. Blocked elements may
be included by specifying zero conductivity values (all elements of tensor) in those locations (also set volume fraction to zero)