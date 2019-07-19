#Calculation name, requires full path
jobid=pwd()*"/../benchmarks/C-Graphene/Graphene_3x3x1_Phonon"

#The choice of primitive or standard cell does not matter since the
#program will use the SPGLIB library to find the primitive cell.
#Mass potential-type-tag fractional-x y z, potential tag should be consistent with LAMMPS format
C = 12.010
# basis = [C 1 0.0 0.166 0;
#          C 1 0.5 0.333 0;
#          C 1 0.5 0.666 0;
#          C 1 0.0 0.833 0 ];
         # C 1 0.5 0.166 0.5;
         # C 1 0.0 0.333 0.5;
         # C 1 0.0 0.666 0.5;
         # C 1 0.5 0.833 0.5];

#Unit cell - Å
# a=2.456;
# b=4.246;
# c=6.696;
# cell = [ a 0.0 0.0;
#          0.0 b 0.0;
#          0.0 0.0 c];

basis = [ C 1 0.333333 0.666666 0.00000;
          C 1 0.666666 0.333333 0.00000];
cell = [2.4543000000 0.0000000000 0.0000000000;
        -1.2271500000 2.1254861485 0.0000000000;
        0.0000000000  0.0000000000 6.6960000000]

#Replicate used on the identified primitive cell
replicate=[2;2;1];

#Turnoff SPGLIB primitive cell finder
spglib = false;

#Potential details, requires full path
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/C.tersoff",
             :potname => :tersoff);

#Phonon Calculation
phonon = false;
#phonononly = true;

#Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =  [0 0 0;
            1//2 1//3 0;
            1//3 1//3 0;
            1//3 1//2 0];



#UTS-8 symbols not yet supported by HDF5
qlabels = [:G;
           :K;
           :M;
           :K]

#Number of points between qpoints
qlinspan = 100;
