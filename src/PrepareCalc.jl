__precompile__()

module PrepareCalc

#External
using LinearAlgebra
using HDF5

#Internal
#using Utilities
using SpglibWrapper

export generate_atomic_structure
export write_unitcellmap_hdf5

@doc """
Return the generated extended cell after finding the primitive cell and creating a basis and unit cell id map.
cell,basis::Array{Real,2} cell and positions of input atomic structure

mass::Array{Real,1} masses of atoms

pottypes::Array{Int,1} the potential type-id, need for energy/force/Hessian routines

replicate::Vector{Int} number of replications in each Cartesian coordinate

TESTED  ✓
""" function generate_atomic_structure(cell::Array{T,2},
                                     basis::Array{T,2},
                                     mass::Array{T,1},
                                     pottypes::Array{Int,1};
                                     replicate::Vector{Int}=[1;1;1],
                                     spglib::Bool=true) where T<:Real
    #Check that basis is scaled coordinates
    if true in (basis .> 1.0)
        error("Basis not in fractional/scaled coordinates.")
    end
    
    
    #First find primitive cell
    mass2pot = Dict(zip(mass,pottypes));

    if spglib == true
        pcell,pbasis,ptypes = find_primitive(cell,basis,mass);
    else
        @info "SPGLIB not used to find primitive cell"
        pcell = cell;
        pbasis = basis;
        #Tuple (mass,basis-id/type)
        ut = unique(mass);
        tid = [findfirst(isequal(u),ut) for u in mass];
        ptypes = [(mass[j],i) for (i,j) in enumerate(tid)];
    end
    
        
    #pcell - Array 3x3
    #pbasis - Array positions(scaled) natoms x 3 
    #ptypes - Tuple(mass,basis-id/type)


    
    spgnum,spgsymb = get_international(pcell,pbasis,ptypes);
    spgstring = String(spgsymb);
    println("------------- Crystal Spacegroup ---------------")
    println("                 #$(spgnum)                     ")
    println("                $(spgstring)                    ")
    println("------------------------------------------------")
    
    npbasis, = size(pbasis);

    
    #Build structure
    nx,ny,nz = replicate[1],replicate[2],replicate[3];

    #Duplicate atoms
    totalcells = nx*ny*nz;
    atoms = Array{Float64,2}(undef,npbasis*totalcells,5);

    #maps unit cell 3D array to an array of atom-id and basis-id tuples
    #example unitcellmap[1,1,1] returns [ (1,1) (1,2) ] 
    unitcellmap = Array{Array{Tuple{Int,Int},1},3}(undef,nx,ny,nz)
    n = 0;
    for x=0:nx-1
         for y=0:ny-1
             for z=0:nz-1
                 extend = [x; y; z];
                 #Array of tuples that are the atom-id and basis-id
                 map = Array{Tuple{Int,Int},1}(undef,npbasis)
                 for b=1:npbasis
                     n = n + 1;
                     atoms[n,1] = ptypes[b][1]; #mass
                     atoms[n,2] = mass2pot[ptypes[b][1]]; #type for potential
                     atoms[n,3:end] = pcell*(pbasis[b,:] + extend);                        
                     map[b] = (n,ptypes[b][2]) #basis-id/type

                 end
                 unitcellmap[x+1,y+1,z+1] = map;
             end
         end
    end
    extendedcell = transpose(pcell) .* replicate;
    return unitcellmap,extendedcell,atoms
end

@doc """
Write the unitcellmap to an HDF5 file.

jobid::String

unitcellmap::AbstractArray
""" function write_unitcellmap_hdf5(jobid::String,
                                    unitcellmap::Array{Array{Tuple{Int,Int},1},3})
    nx,ny,nz = size(unitcellmap)
    nb = length(unitcellmap[1,1,1]);
    nrows= nx*ny*nz*nb;
    ncols = 5; # x y z atom-id basis-type/id
    maparray = zeros(Int64,nrows,ncols);
    n = 0;
    for ix=1:nx
        for iy=1:ny
            for iz=1:nz
                for ib=1:nb
                    n = n+1;
                    maparray[n,1] = ix;
                    maparray[n,2] = iy;
                    maparray[n,3] = iz;
                    # index cell, atom in cell, atom id
                    maparray[n,4] = unitcellmap[ix,iy,iz][ib][1];
                    # index cell, atom in cell, basis type/id
                    maparray[n,5] = unitcellmap[ix,iy,iz][ib][2];
                end #for ib=1:nb
            end
        end
    end #for ix=1:nb

    #Write Map
    wfile = jobid*".hdf5"

    if isfile(wfile)
        println("The file $(wfile) will be overwritten!")
    end
    
    h5open(wfile,"w") do file
        group = g_create(file,"CellMap");
        group["map"] = maparray;
        attrs(group["map"])["Description"] = "Unit cell maping to atom-id and basis type/id.";
        attrs(group["map"])["Columns"] = "Index-x Index-y Index-z Atom-ID Basis-Type/ID"
        attrs(group["map"])["Format"] = "N X 5 Matrix where N is the total number of permutations."
        attrs(group["map"])["Replicate"] = [nx ny nz];
        attrs(group["map"])["NB"] = nb;

    end
end
    
end #GenerateStructure


        

    
