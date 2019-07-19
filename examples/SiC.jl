#Calculation name, requires full path
jobid=pwd()*"/../benchmarks/SiC/SiC_4.32A_2x2x2_Phonon"

#The choice of primitive or standard cell does not matter since the
#program will use the SPGLIB library to find the primitive cell.
#Mass potential-type-tag fractional-x y z, potential tag should be consistent with LAMMPS format
basis = [ 28.086 1 0.0 0 0;
     28.086 1 0.0 0.5 0.5;
     28.086 1 0.5 0.0 0.5;
     28.086 1 0.5 0.5 0.0;
     12.010 2 0.25 0.25 0.25;
     12.010 2 0.25 0.75 0.75;
     12.010 2 0.75 0.25 0.75;
     12.010 2 0.75 0.75 0.25];

#Unit cell - Å
a=4.32;
cell = [ a 0.0 0.0;
         0.0 a 0.0;
         0.0 0.0 a];

#Replicate used on the identified primitive cell
replicate=[2;2;2];

#Potential details, requires full path
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/SiC.tersoff",
             :potname => :tersoff);

#Phonon Calculation
phonon = true;
phonononly = false;

#Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =  [0 0 0;
            1//2 0 1//2;
            3//8 3//8 3//4;
            0 0 0;
            1//2 1//2 1//2] # Γ-X-K-Γ-L

#UTS-8 symbols not yet supported by HDF5
qlabels = [:G;
           :X;
           :K;
           :G;
           :L]

#Number of points between qpoints
qlinspan = 100;
