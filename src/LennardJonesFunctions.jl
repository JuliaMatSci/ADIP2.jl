module LennardJonesFunctions

export getpair

"""
Return Lennard-Jones pairwise interaction:

`` 4\\epsilon \\left(\\left(\\frac{\\sigma}{r}\\right)^{12} - \\left(\\frac{\\sigma}{r}\\right)^{6}\\right ``

r::Array{Float64,1} interaction distance vector i-j

params::Dict{Symbol,T} interaction parameters
""" function getpair(r::Array{N,1},
        params::Dict{Symbol,T}) where {N,T <: Real}
    σ = params[:σ];
    ϵ = params[:ϵ];
    rm = vnorm(r);
    expr = 4ϵ*((σ/rm)^12 - (σ/rm)^6);
    return expr
end #getpair

@doc """
3D Vector norm (L₂ norm/ Euclidean)

vnorm(a)

a::AbstractVector
""" vnorm(a::AbstractVector) = √(a[1]*a[1] +
                              a[2]*a[2] +
                              a[3]*a[3])

end #LennardJonesFunctions
