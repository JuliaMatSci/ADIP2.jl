# ADIP<sup>2</sup>.jl
## Automatic Differentiation of Interatomic Potentials with Phonons in Julia

This is a proof of concept code to see the capability and peformance in using automatic differentiation through the ForwardDiff.jl or ReverseDiff.jl packages. The main use case for this program is to be able to calculate the atomic forces and force constants (i.e. Hessian matrix) for an energy function/routinecall, e.g., lattice sums to get the total energy of a structure. This can then be used to obtain useful lattice properties like phonon dispersion and elastic constants.

The major challenge with this code is the restriction that the function passed to ForwardDiff.jl/ReverseDiff.jl must be unitary (i.e., accept a single argument). This doesn't work well when trying to structure a code to reduce number of out-of-scope calls in a functon for performance considerations. The approach that I've choosen is to pass an array called ``atoms`` which contains the atomic potential-ids/types and spatial coordinates. This however has several issues, 1.) the returned Hessian is contains meaningless zero entries (i.e. due to ids) and 2.) the interatomic potential parameters and cell dimensions cannot be passed or generated within the energy function call, this means they need to be defined globally and done in the code as `` global const ``, 3.) a neighborlist cannot be passed and thus unnessecary looping is performed. It appears as a consequence of these limitations the Hessian calculation is significantly slower than something like a finite-difference approach, and is not practical for large cells. I have still haven't tested, heuristically that is, if defining global variables/hash tables is better than the approach here.

## Dependencies
The standard library LinearAlgebra is required.

There are two major julia packages that are required to run ADIP2.jl:

[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)
[HDF5.jl](https://github.com/JuliaIO/HDF5.jl)

The other required package is a C-API/library called SPGLIB, which provides the routines to determine the primitive unit cell, spacegroup ID, and irreducible reciprocal space mesh. The interface julia module is called SpglibWrapper.jl and within this module
the user needs to change the SPGLIB variable if the dynamical library is not located within the extlib/SPGLIB folder.

[SPGLIB](https://atztogo.github.io/spglib/)

## Supported potentials

1. Tersoff (LAMMPS compatible)
2. Lennard-Jones
3. EAM (Sutton analytical form)

### Adding a potential energy function routine

The process of adding a new potential is relatively straightfoward:

1. Add your potential to the dictionary in Potentials.jl (ex. :lennardjones => "CallLennardJones.jl");
2. Create a file "CallPotentialName.jl" where PotentialName is self-explaning. In this file create a function with the same name as the key you added in step 1., this function should only take single arguement called atoms which is a *n* x 4 AbstractArray  with *n* being the number of atoms. The function should contain the major routine for calculatin the total potential energy of atoms.

## Usage

To run the code from terminal:

`` julia ADIP2.jl path-to-yourscript.jl ```

where `` yourscript.jl `` contains information about the structure and potential. There are three ``REQUIRED`` variables that should be in `` yourscript.jl `` they are:

* ``basis`` is a Nx5 Array with columns: mass id x y z, coordinates need to be in fractional dimensions.
* ``cell`` is a 3x3 Array with : ax ay az  in Angstroms
                                 bx by bz 
				 cx cy cz
* ``replicate`` is the number of replications in each dimensions after the program finds the primitive unit cell.				 
* ``potinfo`` is a dictionary with keyword symbols `` :paramfile `` and ``:potname`` and values that are `` String `` and `` :Symbol `` data types, respectively.

### Phonon option
If your intrested in obtaining the phonon dispersion from the calculated Hessian the following variables shoulds be used:

* ``phonon`` Boolean flag for indicated you want the phonon calculation
* ``phonononly`` in the event that you want the phonon calculation after running ADIP
* ``qpoints`` the irreducible Brillouin zone symmetry points to sample from.
* ``qlabels`` labels for the irreducible symmetry points.
* ``qlinspan`` number of points to sample from between each qpoint.


## Output

The main driver code will output an HDF5 file that contains all the structural, energetic, forces, and Hessian details. If you run the phonon option it will add the results to the same HDF5 file and generate regular text files. The phonon analysis will provide the bandstructure and density of states, the latter is still being polished.

## Benchmarks

![Silicon](benchmarks/Silicon/Si_5.43A_3x3x3_Phonon.pdf| width=300)

![Siliver](benchmarks/Silver/Ag_4.08A_3x3x3_Phonon.png | width=300)

## ISSUES

Currently there seems to be an issue with calculating the energy/force/Hessian of non-orthorhombic cells, this is probably a convention issue.

## Contact

If you have any questions or suggestions to improve this code please contact stefanbringuier at gmail dot com
