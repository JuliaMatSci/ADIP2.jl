#Calculation name
jobid=pwd()*"/../benchmarks/Silicon/Si_5.43A_3x3x3_Phonon"

#Mass potential-type-tag fractional-x y z
basis = [ 28.086 1 0 0 0;
     28.086 1 0.0 0.5 0.5
     28.086 1 0.5 0.0 0.5
     28.086 1 0.5 0.5 0.0
     28.086 1 0.25 0.25 0.25
     28.086 1 0.25 0.75 0.75
     28.086 1 0.75 0.25 0.75
     28.086 1 0.75 0.75 0.25];

#Unit cell - Å
a=5.431;
cell = [ a 0.0 0.0;
         0.0 a 0.0;
         0.0 0.0 a];

#Replicate
replicate=[3;3;3];

#Potential details
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/Si.tersoff.2",
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
