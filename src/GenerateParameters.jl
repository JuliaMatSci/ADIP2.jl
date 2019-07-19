__precompile__()

module GenerateParameters

using Utilities

export genallparams

@doc """
Generate code and interaction parameters based on function call name.

genallparams(potinfo)

potinfo::Dict
""" function genallparams(potinfo::Dict)
    potname = potinfo[:potname];
    paramfile = potinfo[:paramfile];
    if potname == :tersoff
        gentersoffparams(paramfile);
    elseif potname == :lennardjones
        genlennardjonesparams(paramfile);
    elseif potname == :eamsutton
        genparamseamsutton(paramfile);
    end
end

@doc """
TODO

""" function genlennardjonesparams(paramsfile::String)

    potread = parsepotfile(paramsfile);

    #Get number of species                                                                                                                                                                                                                                                                                           
    elset = Set{String}();
    for entry in potread
        e = split(entry)[1];
        push!(elset,e)
    end
    
    #Make map for element to numerical type
    emap = Dict();
    n = 1
    for e in elset
        emap[e] = n;
        n += 1;
    end

    #Generate params dictionary                                                                                                                                                                                                                                                                                      
    params = Dict();
    nlines = length(potread);
    nn = 1;
    flag = false;
    while flag != true
        sline = split(potread[nn])
        e1,e2 = sline[1:2];
        i,j = emap[e1],emap[e2]
        ϵ,σ = parse.(Float64,sline[3:4]);
        pair =  Dict(:σ=>σ,:ϵ=>ϵ);
        if !haskey(params,(i,j))
            params[(i,j)] = pair;
        end
        nn += 1;
        if nn > nlines
            flag = true;
        end
    end
    
    return params
end

@doc """
This parses a potential files for the potential style from 
A.P Sutton, J. Chen Philos. Mag. Lett. 61, 1990. The format is as follows

`` \\rho = \\left(\\frac{a}{r_{ij}}\\right)^{m} ``
`` \\phi = \\left(\\frac{b}{r_{ij}}\\right)^{n} ``

Units are in Å and eV

flag 1
element-1 element-2 pair-exponent pair-term rho-exponent rho-term
... ... ... ... ... ...

TODO: Add cutoff parameter read

""" function genparamseamsutton(paramsfile::String)
    @warn "NOT TESTED."

    potread = parsepotfile(paramsfile);

    #Get EAM flag
    flagid = parse(Float64,potread[1])
    
    #Get number of species                                                                                                                                                                                                                                                                                           
    elset = Set{String}();
    for entry in potread[2:end]
        e = split(entry)[1];
        push!(elset,e)
    end
    
    #Make map for element to numerical type
    emap = Dict();
    n = 1
    for e in elset
        emap[e] = n;
        n += 1;
    end

    #Generate params dictionary                                                                                                                                                                                                                                                                                      
    params = Dict();
    nlines = length(potread);
    nn = 2; #offset from flagid read
    flag = false;
    while flag != true
        sline = split(potread[nn])
        e1,e2 = sline[1:2];
        i,j = emap[e1],emap[e2]
        n,b = parse.(Float64,sline[3:4]);
        m,a = parse.(Float64,sline[5:6]);
        eam =  Dict(:n=>n,:b=>b,:m=>m,:a=>a,
                     :flag=>flagid);
        
        if !haskey(params,(i,j))
            params[(i,j)] = eam;
        end
        nn += 1;
        if nn > nlines
            flag = true;
        end
    end

    return params
end

@doc """
Parse lammps style tersoff file

gentersoffparams(paramsfile)

paramsfile::String
""" function gentersoffparams(paramsfile::String)
    file = open(paramsfile,"r");
    potread = String[];
    for line in readlines(file)
        line = stripcharaddspac(line,'\t');

        if occursin('#',line)
            continue;
        elseif isempty(line)
            continue;
        else
            push!(potread,line)
        end
    end
    close(file);

    #Get number of species
    elset = Set{String}();
    for entry in potread
        e = split(entry)[1];
        push!(elset,e)
    end

    #Make map for element to numerical type
    emap = Dict();
    n = 1
    for e in elset
        emap[e] = n;
        n += 1;
    end


    #Generate params dictionary
    params = Dict();
    nlines = length(potread);
    nn = 1;
    flag = false;
    while flag != true
        sline = split(potread[nn])
        e1,e2,e3 = sline[1:3];
        i,j,k = emap[e1],emap[e2],emap[e3];
        m,γ,λ3 = parse.(Float64,sline[4:6]);
        c,d,h = parse.(Float64,sline[7:9]);
        n,β = parse.(Float64,sline[10:11]);
        λ2,B = parse.(Float64,sline[12:13]);
        R,D = parse.(Float64,sline[14:15]);
        λ1,A = parse.(Float64,sline[16:17]);
        two =  Dict(:A=>A,:B=>B,:λ1=>λ1,:λ2=>λ2,:R=>R,:D=>D,:β=>β,:n=>n);
        three = Dict(:λ3=>λ3,:m=>m,:c=>c,:d=>d,:h=>h,:γ=>γ,:R=>R,:D=>D);
        if !haskey(params,(i,j))
            params[(i,j)] = two;
        end
        params[(i,j,k)] = three;

        nn += 1;
        if nn > nlines
            flag = true;
        end
    end
    return params
end


@doc """
Return the potential file parsed into a string.
    parsepotfile(filename)

filename::String
""" function parsepotfile(filename::String)
    file = open(filename,"r")

    potread = String[];
    for line in readlines(file)
        line = stripcharaddspac(line,'\t');

        if occursin('#',line)
            continue;
        elseif isempty(line)
            continue;
        else
            push!(potread,line)
        end
    end
    close(file);
   return potread
end

end #GenerateParameters
