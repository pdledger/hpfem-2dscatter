# hpfem-2dscatter

This is a simple 2D finite element electromagnetic scattering program that employs a hp-finite element discretisation for the
solution of the curl curl vector wave equation. Discretisations consisting of triangular or quadrilateral meshes are possible
with elements of arbitrary order. Program currently setup for purely triangular meshes for simplicity.

The program has been developed as a MATLAB teaching tool to allow students to play with different types of discretisations
and to experiment with h and p refinements for problems with smooth solutions and for problems with singularities associated
with edges and corners. A perfectly matched layer (PML) approach is used to truncate the otherwise unbounded domain. A number
of benchmark problems are included, and the program computes the radar cross section (per unit length) and contours of the
scattered field. 

The hierarchic H(curl) conforming finite element basis functions are based on those proposed in

S. Zaglmayr High Order Finite Elements for Electromagnetic Field Computation, PhD Thesis, Johannes Kepler University Linz,
Austira 2006 https://www.numa.uni-linz.ac.at/Teaching/PhD/Finished/zaglmayr

The general methodology of the program is described in

P.D. Ledger An hp-adaptive finite element procedure for electromagnetic scattering problems, PhD Thesis, Swansea 
University, UK 2002.

P.D. Ledger K. Morgan, O. Hassan and N.P. Weatherill. Arbitrary order edge elements for electromag- netic scattering problems
using hybrid meshes and a PML. International Journal of Numerical Methods in Engineering, 2002; 55: 339–358.
https://doi.org/10.1002/nme.501

For triangular meshes the code calls the Mesh2d matlab mesh generation tool that has been developed by D. Engwirda

D. Engwirda, Locally-optimal Delaunay-refinement and optimisation-based mesh generation, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, http://hdl.handle.net/2123/13148, 2014.

D. Engwirda, Unstructured mesh methods for the Navier-Stokes equations, Honours Thesis, School of Aerospace, Mechanical and
Mechatronic Engineering, The University of Sydney, 2005.
https://www.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation
