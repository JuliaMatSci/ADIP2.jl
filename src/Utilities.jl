__precompile__()

module Utilities

using HDF5

export printbanner
export formathessian,process_raw_hessian
export stripcharaddspac
export read_and_format_unitcellmap_hdf5

@doc """
Strip characters from regular expression
""" stripcharaddspac(s, r) = replace(s, Regex("[$r]") => " ")

@doc """
Returns a 3N x 3N matrix, where N is the number of atoms and 3 is for the number
of spatial dimensions. This is required when using ForwardDiff/ReverseDiff routines 
since the array atoms passed to the interatomic potential calls includes mass and index
dimensions

hessian,atoms::Array{Real,2}

""" function formathessian(hessian::Array{T,2},
                       atoms::Array{S,2}) where {T <: Real, S <: Real}
    #Reformat Hessian to remove meaning less entries
    n,m = size(atoms);
    frmthess = zeros(n*3,n*3)
    offst = 2
    ioffst = [1 2 3] .+ offst;
    for i=1:n
        joffst = [1 2 3] .+ offst;
        ishift = ioffst .- offst*i
        for j=1:n
            tmphes = hessian[ioffst,joffst]
            jshift = joffst .- offst*j
            frmthess[ishift,jshift] = tmphes;
            joffst = joffst[:].+m;
        end
        ioffst = ioffst[:].+ m;
    end
    return frmthess
end

@doc """
Returns a 3N x 3N matrix, where N is the number of atoms and 3 is for the number
of spatial dimensions. This is required when using ForwardDiff/ReverseDiff routines 
since the array atoms passed to the interatomic potential calls includes mass and index dimensions
as well as not being in the standard xixj xiyj xizj ... format. Here i and j are the atomic indices.

hessian,atoms::Array{Real,2}

""" function process_raw_hessian(hessian::Array{T,2},
                         atoms::Array{S,2}) where {T <: Real, S <: Real}
    #Remove id,type zero entries from hessian
    natoms, = size(atoms);
    offst = natoms + 1;
    slimhess = hessian[offst:end,offst:end]
    
    #Reformat based on traditonal format, see Run.jl run() for comments.
    #input Hessian should be symmetric
    hesscopy = deepcopy(slimhess);
    xx = hesscopy[1:natoms,1:natoms];
    xy = hesscopy[1:natoms,natoms+1:2*natoms];
    xz = hesscopy[1:natoms,2*natoms+1:3*natoms];
    yy = hesscopy[natoms+1:2*natoms,natoms+1:2*natoms]; 
    yz = hesscopy[natoms+1:2*natoms,2*natoms+1:3*natoms];
    zz = hesscopy[2*natoms+1:3*natoms,2*natoms+1:3*natoms];

    #Now loop-over atoms to fill reformated Hessian.
    for i=1:natoms
        ii = 3*(i-1)+1:3*i
        for j=1:natoms
            jj = 3*(j-1)+1:3*j
         
            ij_xyz = zeros(3,3);
            ij_xyz[1,1] = xx[i,j]
            ij_xyz[2,2] = yy[i,j]
            ij_xyz[3,3] = zz[i,j]
            ij_xyz[1,2] = xy[i,j]
            ij_xyz[2,1] = xy[i,j]
            ij_xyz[1,3] = xz[i,j]
            ij_xyz[3,1] = xz[i,j]
            ij_xyz[2,3] = yz[i,j]
            ij_xyz[3,2] = yz[i,j]
            
            slimhess[ii,jj] = ij_xyz;
        end
    end

    return slimhess
end 


@doc """
Read and format the unitcellmap to an HDF5 file. The format will be an Array that is ix X iy X iz with tuple entries.

TESTING: appears to work.
""" function read_and_format_unitcellmap_hdf5(filename::String)
    rfile = h5open(filename,"r");
    data = read(rfile["CellMap"])
    attr = h5readattr(rfile.filename,"CellMap/map")
    nb = attr["NB"];
    nx,ny,nz = attr["Replicate"];
    cellmap = data["map"];

    unitcellmap = Array{Array{Tuple{Int,Int},1},3}(undef,nx,ny,nz)
    n = 0;
    for ix=1:nx
        for iy=1:ny
            for iz=1:nz
                basis = Array{Tuple{Int,Int},1}(undef,nb)
                for ib=1:nb
                    n = n+1;
                    id,btyp = cellmap[n,4:5];
                    basis[ib] = (id,btyp);
                end #for ib=1:nb
                iix,iiy,iiz = cellmap[n,1:3];
                unitcellmap[iix,iiy,iiz] = basis;
            end
        end
    end #for ix=1:nb
    close(rfile);
    return unitcellmap
end

@doc """
Print program banner info.
""" function printbanner()
    println("----------------------------------------------------")
    println("Automatic Differentiation of Interatomic Potentials ")
    println("for Phonons in Julia: ADIP².jl")
    println(" ")
    println("  _____                     _____")
    println(" /     \\       /\\     / \\  /     \\   ")
    println("|       |-    /  \\   /   -|       |   ")
    println(" \\_____/   \\/     \\ /      \\_____/    ")
    println(" ")
    println("Interatomic Potentials:")
    println(" ")
    println("Lennard-Jones ✅")
    println("Tersoff ✅")
    println("EAM (Sutton Analytical) ✅")
    println(" ")
    println("Units:")
    println("Distances in Å")
    println("Energy in eV")
    println("Mass in Kg")
    println(" ")
    println("Limitations: Many, still working through issues.!")
    println(" ")
    println("Author: Stefan Bringuier")
    println("--------------------------------------------------")
end

end
