# monodomain2D
Simple MATLAB code for solution of the Monodomain equation in 2D. Assumes a regular grid, but allows for user to specify heterogeneous diffusion tensors and blocked regions.

Code is commented, so the majority of information is provided there. Use the "create" file as a template to see how to create your own problems to solve. Completely blocked elements
are created simply by specifying zero conductivity values (all elements of tensor) in those locations. New cell models will be added soon, as well as the functionality of specifying
different cell models in different elements.

Note that a version of the monodomain equation allowing for "partial occlusion" by considering the intrinsic average and including volume fraction terms in the equation is what is
solved here. For a standard implementation, simply set volume fraction to one everywhere (and zero for occluded sites). This is also done in the example file. 

Information (mostly the background mathematics) is provided in monodomain2D_FVM.pdf