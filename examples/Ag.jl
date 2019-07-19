#Calculation name
jobid=pwd()*"/../benchmarks/Silver/Ag_4.08A_3x3x3_Phonon"

#Mass potential-type-tag fractional-x y z
Ag = 107.86 
basis = [ Ag 1 0 0 0;
     Ag 1 0.0 0.5 0.5;
     Ag 1 0.5 0.0 0.5;
     Ag 1 0.5 0.5 0.0];

#Unit cell - Å
a=4.09;
cell = [ a 0.0 0.0;
         0.0 a 0.0;
         0.0 0.0 a];

#Replicate
replicate=[3;3;3];

#Potential details
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/Ag.sutton.eam",
             :potname => :eamsutton);

#Phonon Calculation
phonon = true;
phonononly = true;

#Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =  [0 0 0;
            1//2 0 1//2;
            1//2 1//4 3//4;
            1//2 0 1//2;
            3//8 -1//4 3//8;
            0 0 0;
            1//2 1//2 1//2] # Γ-X-W-X-K-Γ-L

#UTS-8 symbols not yet supported by HDF5
qlabels = [:G;
           :X;
           :W;
           :X;
           :K;
           :G;
           :L];

#Number of points between qpoints
qlinspan = 100;
