#Calculation name
jobid=pwd()*"/../benchmarks/Argon/Ar_5.25A_3x3x3_Phonon"

#Mass potential-type-tag fractional-x y z
Ar = 39.948;
basis = [ Ar 1 0 0 0;
          Ar 1 0.0 0.5 0.5;
          Ar 1 0.5 0.0 0.5;
          Ar 1 0.5 0.5 0.0]
    

#Unit cell - Å
a=5.25;
cell = [ a 0.0 0.0;
         0.0 a 0.0;
         0.0 0.0 a];

#Replicate
replicate=[3;3;3];

#Potential details
potinfo=Dict(:paramfile => pwd()*"/../examples/potentialfiles/Ar.lennardjones",
             :potname => :lennardjones);

#Phonon Calculation
phonon = true;
phonononly = false;

#Spacegroup 225 (Fm-3m)
#Brillouin zone sampling (see http://www.cryst.ehu.es/cryst/get_kvec.html)
qpoints =[ 0 0 0;
           1//2 0 1//2;
           1//2 1//2 0;
           0 0 0
           3//8 3//8 0]; #Γ-X-K-Γ-M
            

#UTS-8 symbols not yet supported by HDF5
qlabels = [:G;
           :X;
           :M;
           :G;
           :K];

#Number of points between qpoints
qlinspan = 100;
