#Calculation name, requires full path
jobid=pwd()*"/../benchmarks/GaN/GaN_a3.21A_c5.34A_2x2x2_Phonon"

#The choice of primitive or standard cell does not matter since the
#program will use the SPGLIB library to find the primitive cell.
#Mass potential-type-tag fractional-x y z, potential tag should be consistent with LAMMPS format
Ga = 69.723;
N = 14.006;
basis =[Ga 1 0.333 0.667 1.00;
        Ga 1 0.667 0.333 0.5;
        N 2 0.333 0.667 0.375;
        N 2 0.667 0.333 0.875]

#Unit cell - Å
a1 = [3.216290 0.0 0.0]
a2 = [-1.608145 2.785389 0.0]
a3 = [0.0 0.0 5.239962]
cell = [ a1;
         a2;
         a3];

#Replicate used on the identified primitive cell
replicate=[2;2;2];

#Potential details, requires full path
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/GaN.tersoff",
             :potname => :tersoff);

#Phonon Calculation
phonon = true;
phonononly = true;

#Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =  [0 0 0;
            1//2 0 0;
            1//3 1//3 0;
            0 0 0;
            0 0 1//2;
            1//2 0 1//2;
            1//3 1//3 1//2];

#UTS-8 symbols not yet supported by HDF5
qlabels = [:G;
           :M;
           :K;
           :G;
           :A;
           :L;
           :H]

#Number of points between qpoints
qlinspan = 100;
