using LennardJonesFunctions

@doc """
Main driver function for Lennard-Jones interaction calculation. Does not use
neighbor list and is truncated based on cutoff.

NOTICE: this function requires closure with global constant CELL and
PARAMETERS, the former is instantiated by PrepareCalc.jl and the latter 
by GenerateParameters library call in RunIAP.jl. This is what enables ForwardDiff.jl
calls because at present only unary target function are supported. This has significant
drawbacks and optimal coding structure since its makes use of neighbor list impractical 
given they would have to be declared as globally as  hash maps, although I haven't tested
the performance changes for such approach.

""" function lennardjones(atoms::AbstractArray)
    natoms, = size(atoms);
    types = atoms[:,1];
    pos = atoms[:,2:end];

    #Min. Image Convention
    cellinv = inv(CELL) :: Array{Float64,2};
    spos = pos * cellinv;

    energy = 0.0e0;

    #Complete Loop over atoms
    for i=1:natoms
        for j=1:natoms
            if i == j
                continue
            end

            ityp,jtyp = types[i],types[j];
            ijparams = PARAMETERS[(ityp,jtyp)] :: Dict{Symbol,Float64};

            #rij = pos[j,:]-pos[i,:];
            #Min. image convention for PBC
            srij = spos[j,:] - spos[i,:]
            srij -= round.(srij);
            rij = CELL * srij ;

            vij = getpair(rij,ijparams);

            energy += vij;

            end #j
        end#i
    return energy/2.00e0
end # lennardjones
