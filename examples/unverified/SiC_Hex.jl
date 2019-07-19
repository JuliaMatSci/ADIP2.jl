#Calculation name, requires full path
jobid=pwd()*"/../benchmarks/SiC_a3.09_c15.1_2x2x2_Phonon"

#The choice of primitive or standard cell does not matter since the
#program will use the SPGLIB library to find the primitive cell.
#Mass potential-type-tag fractional-x y z, potential tag should be consistent with LAMMPS format
Si = 28.086;
C = 12.010;
basis = [ Si 1 0.000000 0.000000 0.499723;
Si 1 0.000000 0.000000 0.999723;
Si 1 0.333333 0.666667 0.166509;
Si 1 0.666667 0.333333 0.666509;
Si 1 0.333333 0.666667 0.832986;
Si 1 0.666667 0.333333 0.332986;
C 2 0.000000 0.000000 0.374352 ;
C 2 0.000000 0.000000 0.874352 ;
C 2 0.333333 0.666667 0.041526 ;
C 2 0.666667 0.333333 0.541526 ;
C 2 0.333333 0.666667 0.708004 ;
C 2 0.666667 0.333333 0.208004];

#Unit cell - Å
cell = [ 3.094884 0 0;
         -1.547442 2.680248 0;
          0 0 15.184531];

#Replicate used on the identified primitive cell
replicate=[2;2;2];

#Potential details, requires full path
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/SiC.tersoff",
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
