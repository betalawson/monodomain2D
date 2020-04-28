# monodomain2D

### Code Features
A MATLAB implementation for solution of the [monodomain model](https://en.wikipedia.org/wiki/Monodomain_model), in two dimensions. The code currently uses a regular tetrahedral grid, but allows for the user to specify different diffusion tensors in different elements. Non-diagonal tensors and zero tensors (to represent wholly blocked regions, for example fibrotic obstructions due to materials such as collagen) are also supported. Different cell models may be specified for different nodes, but the current implementation only includes the reduced ten Tusscher *et al.* model (TT3 model). The monodomain model that is solved is actually a modified one that allows for "partial occlusion" by considering the intrinsic average of the typical monodomain equation, leading to a volume fraction term in the equation. The standard monodomain model is recovered simply by using a volume fraction of unity in all unoccupied sites (and zero in occupied sites). The example problem files currently do this.

### Numerical Approach
The code uses a control volume finite element method (CV-FEM), in which a vertex-centred finite volume approach is used, but with interpolants used to better approximate integrals over the control volume, giving an FEM-like quality. The user may select between first and second order timestepping:  
**First Order:** Semi-implicit scheme - Implicit handling of diffusive term, and Rush-Larsen handling of reaction term.  
**Second Order:** Crank-Nicholson handling of diffusive term, Adams-Bashforth and second order generalisation of Rush-Larsen method for reaction term (Perego and Veneziani, 2009).  
The user may also specify to what extent mass lumping should be applied.

### Usage
Code is commented, and the majority of information regarding use of the code is provided via these comments. Basic procedure for usage is to use the "create" file(s) as a template to see how to create your own problems in the structure required by the simulator. Blocked elements may be included by specifying zero conductivity values (all elements of tensor) in those locations (also set volume fraction to zero). Fully occupied sites are also marked with the variable *occ_map*.  
Once a problem to be simulated has been constructed using one of these files, the function *simulateMonodomain* is used to perform a simulation. Users should modify this function to collect the outputs they will need.  
The implementation of the TT3 model provided should be used as a template for the incorporation of other cell models.

### References
Perego, M. & Veneziani, A. (2009). An efficient generalization of the Rush-Larsen method for solving elecro-physiology membrane equations. *Elec. Trans. Num. Anal.* 35: 234-256.  
ten Tusscher, K. H. & Panfilov, A. V. (2006). Cell model for efficient simulation of wave propagation in human ventricular tissue under normal and pathological conditions. *Phys. Med. Biol.* 51: 6141-6156.

### Acknowledgements
This code was created with the support of the ARC Centre of Excellence for Mathematical and Statistical Frontiers   [acems.org.au](http://www.acems.org.au)
