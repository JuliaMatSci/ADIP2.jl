module EAMFunctions

using LinearAlgebra

export getsuttonpair,getsuttonrho,getsuttonembed

@doc """
Returns the pairwise interaction for A.P Sutton, J. Chen, Philos. Mag. Lett. 61, 139, 1990. model

getsuttonpair(r,params)

r::Array{Real,1} distance vector between atom i and j

params:: Dict{Symbol,Real} interaction parameters between i and j

""" function getsuttonpair(r::Array{N,1},
    params::Dict{Symbol,T}) where {N,T<: Real}

    n,b = params[:n],params[:b];
    mr = norm(r);
    pair = (b/mr)^n;
    return pair
end #getsuttonpair

@doc """
Returns the density function for A.P Sutton, J. Chen, Philos. Mag. Lett. 61, 139, 1990. model

getsuttonrho(r,params)

r::Array{Real,1} distance vector between atom i and j

""" function getsuttonrho(r::Array{N,1},
                              params::Dict{Symbol,T}) where {N,T<: Real}

    m,a = params[:m],params[:a];
    mr = norm(r);
    ρ = (a/mr)^m;
    return ρ
end #getsuttonrho

@doc """
Returns the embedding energy for embedding function from 
A.P Sutton, J. Chen, Philos. Mag. Lett. 61, 139, 1990. model

getsuttonembed(ρ)

ρ::Real electron density at a given site

""" function getsuttonembed(ρ::T) where T<:Real
    return -sqrt(ρ)
end

end #EAMFunctions
