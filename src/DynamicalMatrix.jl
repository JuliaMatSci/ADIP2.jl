__precompile__()
module DynamicalMatrix

#External packages
using LinearAlgebra

export getrealspacedynmat
export getrecipspacedynmat

@doc """
Return the mass weighted dynamical matrix in real space.

hessian,mass::Array{Real,2}

""" function getrealspacedynmat(hessian::Array{T,2},
                                mass::Array{S,1}) where {T <: Real, S <: Real}
        # Mass weight normalize
        natoms = length(mass);
        dynmat = zeros(Float64,3*natoms,3*natoms);

        ii = 1;
        for i=1:3*natoms
                jj = 1
                for j=1:3*natoms
                        dynmat[i,j] = 1.00/sqrt(mass[ii]*mass[jj]) * hessian[i,j];
                        if j%3 == 0
                                jj += 1;
                        end
                end
                if i%3 == 0
                        ii += 1;
                end
        end

        #Enforce symmetric nature of dynamical matrix
        dynmat = 0.5e0 * (dynmat+transpose(dynmat));

        #Enforce Acoustic sum rule
        for i=1:3*natoms
            dynmat[i,i] = -sum(dynmat[1:end,i])+dynmat[i,i]
        end
    
        return dynmat
end #getrealspacedynmat

@doc """
Mutable composite datatype for dynamical matrix.

ij::Tuple{Int,Int} indices of atoms i and j.

rij::Vector{Float64} distance vector.

d2::Float64 second derivative (Hessian value).
""" mutable struct HessianEntry
    ij::Tuple{Int,Int}
    rij::Vector{Float64};
    d2::Float64;
end

@doc """
Return the mass weighted dynamical matrix in real space.

hessian,mass::Array{Real,2}

""" function getrealspacedynmat(hessian::Array{T,2},
                                cell::Array{T,2},
                                positions::Array{T,2},
                                types::Array{Int,1},
                                mass::Array{S,1}) where {T <: Real, S <: Real}

        natoms = length(mass);
        dynmat = Array{HessianEntry,2}(undef,3*natoms,3*natoms);

        #Min. Image Convention
        cellinv = inv(cell) :: Array{Float64,2};
        spos = positions * cellinv;

        ii = 1;
        for i=1:3*natoms
                jj = 1;
            for j=1:3*natoms
                ij = (ii,jj);
                srij = spos[jj,:] - spos[ii,:];
                srij -= round.(srij);
                rij = cell * srij ;
                massweight_d2drij = 1.00/sqrt(mass[ii]*mass[jj]) * hessian[i,j];
                dynmat[i,j] = HessianEntry(ij,rij,massweight_d2drij); 
                        if j%3 == 0
                                jj += 1;
                        end
                end
                if i%3 == 0
                        ii += 1;
                end
        end
        
        #Enforce symmetric nature of dynamical matrix
        #Lazy loop
        for i=1:3*natoms
            for j=1:3*natoms
                dij = dynmat[i,j].d2;
                dji = dynmat[j,i].d2;
                dynmat[i,j].d2 = 0.5e0 * (dij+dji);
                dynmat[j,i].d2 = 0.5e0 * (dij+dji);
            end
        end
    
       #Enforce Acoustic sum rule
        for i=1:3*natoms
        offsum = 0.00e0
        for j=1:3*natoms
           offsum  += dynmat[i,j].d2
        end
        dynmat[i,i].d2 = dynmat[i,i].d2 - offsum
       end
    return dynmat
end #getrealspacedynmat

@doc """
Return the reciprocal space representation of dynamical matrix

realdynmat::AbstractArray <: Array{HessianEntry,2} mass weight dynamical matrix

unitcellmap::AbstractArray <: Array{Array{Tuple{Int,Int},1},1} atomic mapping
to unitcell id and basis type.

qmesh:: Array{Real,2} reciprocal space mesh sampling

NOTES: 
1. This is by no means an optimal way to construct the reciprocal space representation
of the dynamical matrix, however, its what made sense to me. I WELCOME SUGGESTIONS FOR IMPROVEMENTS.

2. This requires a reciprocal space mesh based on the primitive cell and basis
""" function getrecipspacedynmat(realdynmat::AbstractArray,
                                 unitcellmap::AbstractArray,
                                 qmesh::Array{T,2}) where {T <: Real}

    nentries, = size(realdynmat);
    natoms = convert(Int,nentries//3);
    nx,ny,nz = size(unitcellmap);
    
    # First restructure the dynamical matrix in terms of a dictonary/set
    # for each atom. This means that we are slicing up the dynamical matrix
    # such that a the atom-id will return an Array{HessianEntry,2}(nentries,3)
    drealdynmat = Dict();

    for i=1:natoms
        ii = 3*(i-1)+1:3*i
        drealdynmat[i] = realdynmat[1:end,ii];
    end

    # Now extract the atom-ids and basis-types/ids for the fist unit cell
    # which is an Array{Tuple{Int,Int},1} where each Tuple is atom-id,atom-basis
    firstprimbasis = unitcellmap[1,1,1]; 
    
    # Incorporate values of firt unit cell for each q-vector
    nbasis = length(firstprimbasis);
    nq, = size(qmesh);
    qdynmat = zeros(ComplexF64,3*nbasis,3*nbasis,nq);

    btypposmap = Dict();
    #Wave vector loop
    for qi=1:nq
        q = qmesh[qi,:];
        dq = zeros(ComplexF64,3*nbasis,3*nbasis)
        
        #Loop over first unit cell primitive basis atoms
        for p in firstprimbasis
            pid,pbtyp = p[1],p[2]  
            ipid = 3*(pid-1)+1:3*pid;
            for pp in firstprimbasis
                ppid,ppbtyp = pp[1],pp[2];
                ippid = 3*(ppid-1)+1:3*ppid;

                #Extract the appropriate dynamical matrix entries
                dynmat_entry = drealdynmat[pid][ippid,1:end] #This should be Array{HessianEntry,2}(3,3)
                @assert size(dynmat_entry) == (3,3)
                #Get ids and check tye match
                ab=dynmat_entry[1,1].ij;
                
                #@assert ab[1] == pid && ab[2] == ppid
                
                #Get local mass weigh Hessian
                tmpdyn = zeros(Float64,3,3);
                for a=1:3
                    for b=1:3
                        tmpdyn[a,b] = dynmat_entry[a,b].d2;
                    end
                end

                #Get distance position between ipid and ippid
                #doesn't matter which point we access
                r = dynmat_entry[1,1].rij;

                #Create map of basis type to position useful for next major set of loops
                btypposmap[(pbtyp,ppbtyp)] = r;
                
                #Construct transform
                dq[ipid,ippid] = tmpdyn[:,:] * exp(-im*dot(r,q));
            end
        end
        qdynmat[:,:,qi] = dq[:,:];
    end

    #Now accumulate over all unit cells for each wavevector
     for qi=1:nq
         q = qmesh[qi,:];
         #Start with primitive reference unit cell.
         for p in firstprimbasis
             pid,pbtyp = p[1],p[2];
             ipid = 3*(pid-1)+1:3*pid;
             #Loop over all other unit cells
             for ix=1:nx
                 for iy=1:ny
                     for iz=1:nz
                         #Skip primitive reference unit cel
                         if ix == 1 && iy == 1 && iz == 1
                             continue
                         end
                         #Get atom unit cell mapping
                         uat = unitcellmap[ix,iy,iz]; #Array{Tuple{Int,Int}}

                         #Get reciprocal space dynmat for atom pairs and
                         #for the correct basis site in primitive unit cell.
                         for a in uat
                             id,btyp = a[1],a[2]
                             iid = 3*(id-1)+1:3*id
                             iibtyp = 3*(btyp-1)+1:3*btyp;
                             dynmat_entry = drealdynmat[pid][iid,1:end];
             
                             #Get local mass weighted Hessian                                                                                                                                                      
                             tmpdyn = zeros(Float64,3,3);
                             for a=1:3
                                 for b=1:3
                                     tmpdyn[a,b] = dynmat_entry[a,b].d2;
                                 end
                             end
                             #Get distance position between ipid and ippid                                                                                                                                      
         	             #doesn't matter which point we access                                                                                                                                              
                             #Usef basis-type to position map to shift back to primitive unit cell
                             r = dynmat_entry[1,1].rij;

                             #Construct transform                                                                                                                                                               
                             qdynmat[ipid,iibtyp,qi] += tmpdyn[:,:] * exp(-im*dot(r,q));
                         end #for a in uat
                     end
                 end
             end #for ix=1:nx
         end # for p in firstprimbasis
     end #for qi=1:nq

    return qdynmat
end


end #DynamicalMatrix

    
        
