module SpglibWrapper

#External
using LinearAlgebra


###USER CHANGE if SPGLIB is in different path
const SPGLIB = joinpath(dirname(@__FILE__),"../extlib/SPGLIB/_build/libsymspg.so")
###

export find_primitive
export get_international
export get_ir_reciprocal_mesh

@doc """
Spglib wrapper to return the primitive unit cell for lattice R⃗ with atomic positions r⃗ of type t.

R⃗,r⃗::Array{Real,2}

t::AbstractArray

symprec::Real

REQUIREMENTS: This function realies on a call to the spg_standardize_cell SPGLIB library call. The
library can be found at https://atztogo.github.io/spglib/. A global const spglib is defined in this
module which contains the path to the library.

NOTES: Assumes input cell and basis is in the format:
cell =  ax ay az
        bx by bz
        cx cy cz
basis x₁ y₁ z₁
      x₂ y₂ z₃ 
      .....

""" function find_primitive(R⃗::Array{T,2},r⃗::Array{T,2},mass::AbstractArray;
                        symprec::T=1.0e-5) where {T <: Real}

    #Spglib appears to want cell basis vectors column major.
    #Spglib appears to want position x atoms-id format
    cR⃗ = convert(Array{Cdouble},transpose(deepcopy(R⃗))); # CHECK: should this be transposed
    cr⃗ = convert(Array{Cdouble},transpose(deepcopy(r⃗)));
    ut = unique(mass)

    #The type array is given as Masses
    mt = Dict();
    for i=1:length(ut)
        mt[i]=ut[i]
    end
    
    tᵢ = convert(Array{Cint},[findfirst(isequal(u),ut) for u in mass]);
    natoms = length(tᵢ);

    natomprim = ccall(
    (:spg_standardize_cell,SPGLIB),
    Cint,
        (Ptr{Cdouble},Ptr{Cdouble},
         Ptr{Cint},Cint,
         Cint,Cint,Cdouble),
         cR⃗,cr⃗,tᵢ,
         natoms,1,0,symprec);
    
    if natomprim == 0
        error("Could not determine primitive cell!")
    end
    R⃗′ = convert(Array{Float64},transpose(cR⃗));
    r⃗′ = convert(Array{Float64},transpose(cr⃗)[1:natomprim,1:end]);
    #Return a array with tuple entries of basis type and mass
    ut′ = [(mt[j],i) for (i,j) in enumerate(tᵢ[1:natomprim])];
    
    return (R⃗′,r⃗′,ut′)
end

@doc """
Spglib wrapper to get the spacegroup (international representation) for unit cell with lattice
R⃗ and positions r⃗ of type t.

R⃗,r⃗::Array{Real,2}

t::AbstractArray

symprec::Real

REQUIREMENTS: This function realies on a call to the spg_standardize_cell SPGLIB library call. The
library can be found at https://atztogo.github.io/spglib/. A global const spglib is defined in this
module which contains the path to the library.

NOTES: Assumes input cell and basis is in the format:
cell =  ax ay az
        bx by bz
        cx cy cz
basis x₁ y₁ z₁
      x₂ y₂ z₃ 
      .....

""" function get_international(R⃗::Array{T,2}, r⃗::Array{T,2}, mass::AbstractArray;
                               symprec::T=1.0e-4) where {T <: Real}
    
    spgsymbol = Array{UInt8}(undef,11);
    cR⃗ = convert(Array{Cdouble},transpose(R⃗)); # CHECK: should this be transposed

    #Assumes position array is in  position by atoms-id format

    cr⃗ = convert(Array{Cdouble},transpose(r⃗));
    ut = unique(mass)

    #The type array is given as Masses
    mt = Dict();
    for i=1:length(ut)
        mt[i]=ut[i]
    end
    
    tᵢ = convert(Array{Cint},[findfirst(isequal(u),ut) for u in mass]);
    natoms = length(tᵢ);

    spgid = ccall(
    (:spg_get_international,SPGLIB),
    Cint,
        (Ptr{Cuchar},Ptr{Cdouble},
         Ptr{Cdouble},Ptr{Cint},
         Cint,Cdouble),
         spgsymbol,cR⃗,cr⃗,tᵢ,
         natoms,symprec);
    
    if spgid == 0
        error("Could not determine spacegroup!")
    end

    return (spgid,spgsymbol)
end
                 
@doc """
Not Tested!
Spglib wrapper for returning the reciprocal q-mesh using.
""" function get_ir_reciprocal_mesh(mesh::AbstractVector,
                                lattice::AbstractArray,
                                positions::AbstractArray,
                                types::AbstractVector;
                                symprec::Real=1e-6,
                                center::Array{T} = [0.0 0.0 0.0]) where {T <: Real} 

    mapping = zeros(Cint,prod(mesh));
    grid = zeros(Cint,prod(mesh),3);

    cmesh = convert(Array{Cint},mesh);
    clattice = convert(Array{Cdouble},transpose(lattice));
    cpositions = convert(Array{Cdouble},transpose(positions));
    ccenter = convert(Array{Cdouble},center);
    ctimerev = 0;
    unique_types = unique(types)
    type_indices = Cint[findfirst(isequal(u),unique_types) for u in types]
    natoms = length(type_indices);

    numops = ccall(
    (:spg_get_ir_reciprocal_mesh,SPGLIB),
    Cint,
        (Ptr{Cint},Ptr{Cint},
         Ptr{Cint},Ptr{Cint},
         Cint,Ptr{Cdouble},
         Ptr{Cdouble},Ptr{Cint},
         Cint,Cdouble),
        grid,mapping,
        cmesh,ccenter,
        ctimerev,clattice,
        cpositions,type_indices,
        natoms,symprec)
    
    if numops == 0
    error("Could not determine irreducible grid and mesh")
    end
    #Get unique q-points in mesh. Need to shift for index purposes
    uslice = unique(convert(Array{Int},mapping)).+1
    qgrid = convert(Array{Int},grid);
    qgrid = qgrid[uslice,:]./ [mesh[1] mesh[2] mesh[3]]; 
    return qgrid
end


end #SpglibWrapper.jl
