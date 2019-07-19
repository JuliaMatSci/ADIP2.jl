module TersoffFunctions

using LinearAlgebra

export getcutoff
export getattractive,getrepulsive,getbondorder
export getthreebody,getpairmanybody
export getzeta


@doc """
Returns evaluated cutoff function:

``  0.5 - 0.5\\sin\\left(\\frac{\\pi}{2}\\frac{r-R}{D}\\right) ``

r::Float64 interaction distance

params::Dict{Symbol,Real} interaction parameters
""" function getcutoff(r::Array{N,1},
    params::Dict{Symbol,T}) where {N,T<: Real}

    R,D = params[:R],params[:D];
    mRD = R-D;
    pRD = R+D;
    mr = norm(r);

    if mr < mRD
        return 1
    elseif mr > mRD && mr < pRD
        return 0.5-0.5*sin((pi/2)*(mr-R)/D)
    else
        return 0
    end
end #getcutoff

@doc """
Returns the i-j bondorder term used to modify the attractive potential:

`` \\left(1 + \\beta^{n} \\zeta^{n} \\right)^{\\frac{-1}{2n}} ``

ζ::Float64 three-body accumaltive interaction term

params::Dict{Symbol,Real} interaction parameters
""" function getbondorder(ζ::N,
        params::Dict{Symbol,T}) where {N,T <: Real}

    β,n = params[:β],params[:n];

    expr = (1 + β^n * ζ^n)^(-1/(2*n));
    return expr
end #getbondorder

@doc """
Returns the function product 3-body interactions.

getzeta(r1,r2,params)

r1,r2::Array{Float64,1} interaction distance vectors for i-j and i-k

params::Dict{Symbol,Real} interaction parameters
""" function getzeta(r1::Array{N,1},r2::Array{N,1},
        params::Dict{Symbol,T}) where {N,T <: Real}

    cutfr2 = getcutoff(r2,params);
    gtheta = getthreebody(r1,r2,params);
    expr1r2 = getpairmanybody(r1,r2,params);

    expr = cutfr2*gtheta*expr1r2
    return expr
end #getzeta

"""
Return the i-j-k threebody interaction:

`` \\gamma \\left( 1 + \\frac{c^2}{d^2} - \frac{c^2}{d^2 + \\left(\\cos \\theta - h\\right)^2}\\right) ``

getthreebody(r1,r2,params)

r1,r2::Array{Float64,1} interaction distance vectors i-j and i-k

params::Dict{Symbol,Real} interaction parameters
""" function getthreebody(r1::Array{N,1},r2::Array{N,1},
             params::Dict{Symbol,T}) where {N,T <: Real}

    c,d,h = params[:c],params[:d],params[:h];
    γ = params[:γ];
    costheta = dot(r1,r2)/(norm(r1)*norm(r2));
    term1 = (c/d)^2;
    term2 = c^2 / ( d^2 + (costheta-h)^2);

    expr = γ*(1 + term1 - term2);
    return expr
end #getthreebody

"""
Return the i-j and i-k 3-body pairwise interaction:

`` \\exp^{\\lambda_{3}^{m} \\left(r_{1}-r_{2}\\right)} ``

r1,r2::Array{Float64,1} interaction distance vectors i-j and i-k

params::Dict{Symbol,Real} interaction parameters
""" function getpairmanybody(r1::Array{N,1},r2::Array{N,1},
             params::Dict{Symbol,T}) where {N,T <: Real}

    λ3,m = params[:λ3],params[:m];
    mr1 = norm(r1);
    mr2 = norm(r2);

    expr = exp(λ3^m * (mr1-mr2)^m);
    return expr
end #getpairmanybody

"""
Return Morse based i-j pairwise attractive interaction:

``` -B \\exp\\left(-\\lambda_{2} r\\right) ```

r1::Array{Float64,1} interaction distance vector i-j

params::Dict{Symbol,Real} interation parameters
""" function getattractive(r::Array{N,1},
        params::Dict{Symbol,T}) where {N,T <: Real}
    B = params[:B];
    λ2 = params[:λ2];
    mr = norm(r);

    expr = -B*exp(-λ2*mr);
    return expr
end #getattractive

"""
Return the Mose based i-j pairwise repulsive interaction.

`` A\\exp\\left(-\\lambda_{1} r \\right) ``

r::Array{Float64,1} interaction distance vector i-j

params::Dict{Symbol,Real}
""" function getrepulsive(r::Array{N,1},
        params::Dict{Symbol,T}) where {N,T <: Real}
    A = params[:A];
    λ1 = params[:λ1];
    mr = norm(r);

    expr = A*exp(-λ1*mr);
    return expr
end #getrepulsive

end #TersoffFunctions
