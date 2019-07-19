
a,c = 5.0563,7.3739
Nx,Ny,Nz = 1,1,1;
#Mass type x y z
basis = [28.086 1 0.000000 0.000000 0.500000;
         28.086 1 0.000000 0.500000 0.750000;
         28.086 1 0.500000 0.500000 0.000000 ;
         28.086 1 0.500000 0.000000 0.250000 ;
         16.00 2 0.087966 0.750000 0.625000 ;
         16.0  2 0.750000 0.912034 0.375000 ;
         16.0  2 0.912034 0.250000 0.625000 ;
         16.0  2 0.250000 0.087966 0.375000 ;
         16.0  2 0.250000 0.412034 0.875000 ;
         16.0  2 0.750000 0.587966 0.875000 ;
         16.0  2 0.587966 0.250000 0.125000 ;
         16.0  2 0.412034 0.750000 0.125000];
nbatoms, = size(basis)
cell = [ a*Nx 0.0 0.0;
                0.0 a*Ny 0.0;
                0.0 0.0 c*Nz];

#Duplicate atoms
atoms = Array{Float64,2}(undef,length(basis[:,1])*Nx*Ny*Nz,5);
let
    n = 0;
    for x=0:Nx-1
         for y=0:Ny-1
             for z=0:Nz-1
                 for b=1:nbatoms
                     n = n + 1;
                     atoms[n,1] = basis[b,1];
                     atoms[n,2] = basis[b,2];
                     atoms[n,3] = a*(x+basis[b,3]);
                     atoms[n,4] = a*(y+basis[b,4]);
                     atoms[n,5] = a*(z+basis[b,5]);
                 end
             end
         end
     end
 end

potinfo=Dict(:paramfile => "scripts/SiO.tersoff", :potname => :tersoff);
