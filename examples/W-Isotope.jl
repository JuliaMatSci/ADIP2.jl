#Calculation name
jobid=pwd()*"/../benchmarks/Tungsten/WIsotope_3.165A_3x3x3_Phonon"

#Mass potential-type-tag fractional-x y z
W182 = 181.9482042
W186 = 185.954364 
basis = [ W182 1 0 0 0;
          W186 1 0.5 0.5 0.5]; #Slight offset position to get optical bands (W-W dimer)

#Unit cell - Å
a=3.165;
cell = [ a 0.0 0.0;
         0.0 a 0.0;
         0.0 0.0 a];

#Replicate
replicate=[3;3;3];

#Potential details
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/W.Juslin",
             :potname => :tersoff);

#Phonon Calculation
phonon = true;
phonononly = true;

#  Spacegroup #229 Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =  [1//4 1//4 1//4;
            0 0 1//2;
            0 0 0;
            1//2 -1//2 1//2;
            1//4 1//4 1//4;
            0 0 0] # P-N-Γ-H-P-Γ

#UTS-8 symbols not yet supported by HDF5
qlabels = [:P;
           :N;
           :G;
           :H;
           :P;
           :G];

#Number of points between qpoints
qlinspan = 100;
