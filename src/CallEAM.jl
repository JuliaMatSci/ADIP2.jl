using EAMFunctions

@doc """
Main driver function for EAM interaction calculation from
A.P Sutton, J. Chen, Philos. Mag. Lett. 61, 139, 1990. model. Does not use
neighbor list and is truncated based on cutoff.

NOTICE: this function requires closure with global constant CELL and
PARAMETERS, the former is instantiated by PrepareCalc.jl and the latter 
by GenerateParameters library call in RunIAP.jl. This is what enables ForwardDiff.jl
calls because at present only unary target function are supported. This has significant
drawbacks and optimal coding structure since its makes use of neighbor list impractical 
given they would have to be declared as globally as  hash maps, although I haven't tested
the performance changes for such approach.

""" function eamsutton(atoms::AbstractArray)
    natoms, = size(atoms);
    types = atoms[:,1];
    pos = atoms[:,2:end];

    #Min. Image Convention
    cellinv = inv(CELL) :: Array{Float64,2};
    spos = pos * cellinv;


    #Build-up density on atom i and
    # get pairwise energy first
    rho =zeros(eltype(atoms),natoms) #ForwardDiff will complain if you specify Float64
    pairenergy = 0.00e0;
    for i=1:natoms
        ityp = types[i];
        for j=1:natoms
            if i == j
               continue
            end
            sposi = spos[i,:];
            sposj = spos[j,:];            

            #rij = pos[j,:]-pos[i,:];
            #Min. image convention for PBC
            srij = sposj - sposi;
            srij -= round.(srij);
            rij = CELL * srij ;

            jtyp = types[j];
            ijparams = PARAMETERS[(ityp,jtyp)] :: Dict{Symbol,Float64};

            #TODO: Add cutoff check?
            # if norm(rij) >= ijparms[:R]+ijparms[:D]
            #     continue   
            # end
            rho[i] += getsuttonrho(rij,ijparams); 
            pairenergy += getsuttonpair(rij,ijparams);
        end #j
    end#i
    pairenergy /= 2.00e0;

    #Now get embedding energy on atom i
    embedenergy = 0.00e0;
    for i=1:natoms
        embedenergy += getsuttonembed(rho[i]);
    end

    # E = F(ρ) + 1/2 ϕ(r)
    energy = embedenergy + pairenergy;
    return energy
        
end 
