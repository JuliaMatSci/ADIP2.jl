#Calculation name, requires full path
jobid=pwd()*"/../benchmarks/GaN_ZincBlende_4.49A_3x3x3_Phonon"

#The choice of primitive or standard cell does not matter since the
#program will use the SPGLIB library to find the primitive cell.
#Mass potential-type-tag fractional-x y z, potential tag should be consistent with LAMMPS format
Ga = 69.723;
N = 14.006;
basis = [ Ga 1 0.0 0 0;
     Ga 1 0.0 0.5 0.5
     Ga 1 0.5 0.0 0.5
     Ga 1 0.5 0.5 0.0
     N 2 0.25 0.25 0.25
     N 2 0.25 0.75 0.75
     N 2 0.75 0.25 0.75
     N 2 0.75 0.75 0.25];

#Unit cell - Å
a=4.4985;
cell = [ a 0.0 0.0;
         0.0 a 0.0;
         0.0 0.0 a];

#Replicate used on the identified primitive cell
replicate=[3;3;3];

#Potential details, requires full path
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/GaN.tersoff",
             :potname => :tersoff);

#Phonon Calculation
phonon = true;
phonononly = false;

#Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =  [1//2 0 1//2;
            0 0 0;
            1//2 1//2 1//2;
            1//2 0 1//2;
            1//2 1//4 3//4;
            1//2 1//2 1//2] # X-Γ-L-X-W-L

#UTS-8 symbols not yet supported by HDF5
qlabels = [:X;
           :G;
           :L;
           :X;
           :W;
           :L]

#Number of points between qpoints
qlinspan = 100;
