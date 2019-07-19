using TersoffFunctions

@doc """
Main driver function for Tersoff interaction calculation. Does not use
neighbor list and is truncated based on cutoff.

NOTICE: this function requires closure with global constant CELL and
PARAMETERS, the former is instantiated by PrepareCalc.jl and the latter 
by GenerateParameters library call in RunIAP.jl. This is what enables ForwardDiff.jl
calls because at present only unary target function are supported. This has significant
drawbacks and optimal coding structure since its makes use of neighbor list impractical 
given they would have to be declared as globally as  hash maps, although I haven't tested
the performance changes for such approach.

""" function tersoff(atoms::AbstractArray)
    natoms, = size(atoms);
    #mass = atoms[:,1];
    types = atoms[:,1];
    pos = atoms[:,2:end];

    #Min. Image Convention
    cellinv = inv(CELL) :: Array{Float64,2};
    spos = pos * cellinv;

    energy = 0.0e0;

    #Complete Loop over atoms
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
            ijparms = PARAMETERS[(ityp,jtyp)] :: Dict{Symbol,Float64};

            if norm(rij) >= ijparms[:R]+ijparms[:D]
                continue   
            end

            # Eq. $f_c\left(V_r+b_{ij}V_a\right)$
            fcij = getcutoff(rij,ijparms)
            uaij = getattractive(rij,ijparms)
            urij = getrepulsive(rij,ijparms)

            ζij = 0.00e0;
            for k=1:natoms
                if j == k || i == k
                    ζij += 0.00e0;
                else                
                    sposk = spos[k,:];
            
                    #rik = pos[k,:]-pos[i,:];
                    #Min image convention for PBC
                    srik = sposk - sposi;
                    srik -= round.(srik);
                    rik = CELL * srik;
                    
                    ktyp = types[k];
                    ijkparms = PARAMETERS[(ityp,jtyp,ktyp)] :: Dict{Symbol,Float64};

                    if norm(rik) >= ijkparms[:R]+ijkparms[:D]
                        continue
                    end

                    ζij += getzeta(rij,rik,ijkparms);
                end
            end #k

            bij = getbondorder(ζij,ijparms);

            #Tada!!
            vij = fcij*(urij+bij*uaij);

            energy += vij;

            end #j
        end#i
    return energy/2.00e0
end
