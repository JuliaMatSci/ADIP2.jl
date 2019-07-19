__precompile__()
module ReciprocalSpace

#using Utilities
using SpglibWrapper

export construct_reciprocal_space
export construct_ir_mesh

@doc """
Function to generate primitive cell/basis and q-sampling based on provided
high-symmetry q-points (e.g. \\Gamma -> X ) in fractional coordinates. 

R⃗::Array{Real,2} calculation cell in Å

r⃗::Array{Real,2} atomic positions in fractional coordinates

m::Array{Real,1} masses or types of atoms

q⃗::Array{Number,1} fractional q-points in reciprocal space.

REQUIREMENTS:
This function needs access to the function ```find_primitive``` in Utilities.jl which
in turn calls upon the SPGLIB library that is used to determine the primitive cell. For
the SPGLIB library please see details at https://atztogo.github.io/spglib/.

""" function construct_reciprocal_space(R⃗::Array{T,2},
                                        r⃗::Array{T,2},
                                        m::Array{T,1},
                                        q⃗::Array{S,2};
                                        n::Int=100) where {T <: Real, S <: Number}
    #Need to ENFORCE that r⃗ is in scaled coordinates
    if maximum(r⃗) > 1.0e-16
        r⃗ = r⃗ * inv(R⃗)
    end
    R⃗′,r⃗′,t′ = find_primitive(R⃗,r⃗,m);
    q⃗_mesh,q⃗_norm = get_kspace_vectors(q⃗,R⃗′,n=n)
    return (r⃗′,t′,q⃗_mesh,q⃗_norm)
end

@doc """

""" function construct_ir_mesh(R⃗::Array{T,2},
                               r⃗::Array{T,2},
                               m::Array{T,1},
                               qmesh::Array{Int,1}) where {T<:Real}
    R⃗′,r⃗′,t′ = find_primitive(R⃗,r⃗,m);
    return get_ir_reciprocal_mesh(qmesh,R⃗′,r⃗′,m)
end

    
@doc """
Return k-space sampling of vectors q⃗ in primitive cell R⃗ with n number of points.
q⃗,R⃗ :: Array{Real,2}
n :: Int

NOTES: q⃗ should be fractional coordinates.

REQUIREMENTS:
This function needs access to the function ```find_primitive``` in Utilities.jl which
in turn calls upon the SPGLIB library that is used to determine the primitive cell. For
the SPGLIB library please see details at https://atztogo.github.io/spglib/.

""" function get_kspace_vectors(q⃗::Array{S,2},
                                R⃗::Array{T,2};
                                n::Int=100) where {S <: Number, T <: Real}
    num_q⃗, = size(q⃗);
    G⃗ = get_reciprocal_lattice(R⃗);

    q⃗ᵣ= zeros(Float64,num_q⃗,3);
    for i=1:num_q⃗
        q⃗ᵣ[i,:] = G⃗ * q⃗[i,:];
    end

    #Generate q-point mesh
    num_q⃗_points = n * (num_q⃗-1);
    q⃗_mesh = zeros(Float64,num_q⃗_points,3);
    q⃗_norm = zeros(Float64,num_q⃗);
    # Julia is column-major indexing, so this can
    # be slow for a large number of points.
    for i=1:num_q⃗-1
        q⃗_norm[i+1] = L₂(q⃗ᵣ[i+1,:]-q⃗ᵣ[i,:])
        ib = (i-1) * n + 1 ;
        ie = i * n ; 
            for j=1:3
                q⃗_path = range(q⃗ᵣ[i,j],stop=q⃗ᵣ[i+1,j],length=n);
                q⃗_mesh[ib:ie,j] = collect(q⃗_path);
        end
    end

    return q⃗_mesh,cumsum(q⃗_norm)
end

@doc """
3D Vector norm (L₂ norm/ Euclidean)
a::AbstractVector
""" L₂(a::AbstractVector) = √(a[1]*a[1] +
                              a[2]*a[2] +
                              a[3]*a[3])

@doc """
3D Vector cross product a x b
a,b::AbstractVector
""" (×)(a::AbstractVector,b::AbstractVector) = [a[2]*b[3] - a[3]*b[2];
                                                a[3]*b[1] - a[1]*b[3];
                                                a[1]*b[2] - a[2]*b[1]];
@doc """
3D Vector dot(\\cdot) product a . b
a,b::AbstractVector
""" (⋅)(a::AbstractVector,b::AbstractVector) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3];


@doc """
Matrix transpose
a::AbstractMatrix

TODO: make this an infix operator \\^T<tab>

""" bruteT(a::AbstractMatrix) = begin
    n,m = size(a);
    b = zeros(eltype(a),m,n);
    for i=1:n
        for j=1:m
            b[j,i] = a[i,j];
        end
    end
    return b
end

    
@doc """
Return reciprocal space lattice.
""" function get_reciprocal_lattice(R⃗::Array{Float64,2})
    (r₁,r₂,r₃) = R⃗[1,:],R⃗[2,:],R⃗[3,:];
    Ω = r₁⋅(r₂ × r₃)  ;
    b₁ = 2π *  r₂ × r₃ ;
    b₂ = 2π *  r₃ × r₁ ;
    b₃ = 2π *  r₁ × r₂ ;
    G⃗ = [ b₁  b₂  b₃ ] ./ Ω ;
    return bruteT(G⃗)
end




end #ReciprocalSpace
