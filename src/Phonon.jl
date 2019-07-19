__precompile__()
module Phonon

#External
using DelimitedFiles
using LinearAlgebra
using HDF5

#Internal
using Utilities 
using DynamicalMatrix
using ReciprocalSpace
using SpglibWrapper

export phonons

const THz = 15.63408085e0 #Atomic units to THz (THz/AMU)

@doc """ 
Claculate the phonon dispersion from an ADIP.jl calculation. Saves to HDF5 and
returns wavenumber magnitude line and frequency bands.

cell,hessian,positions::Array{Real,2} calculation details

masses::Array{Real,1} 

types::Array{Int,1} atomic potential types

unitcellmap::AbstractArray atomic map for associated unitcell id and basis type

qlinspan::Int number of points between reciprocal space q-space points

qpoints::Array{Number,2} reciprocal (q-space) points, typically high-symmetry.

STEPS PERFORMED BY THIS FUNCTION:
    1.) Calculate the reciprocal mesh for calculation cell/positions
    2.) Generate reciprocal space distance line.
    3.) Get the mass weighted dynamical matrix in real space coordinates.
    4.) Get the reciprocal space representation from step 3 output.
    5.) Calculate the eigenvalues/freq. for each q-point in q-mesh
    6.) Save to HDF5 file and return 

""" function phonons(cell::Array{T,2},
                              hessian::Array{T,2},
                              positions::Array{T,2},
                              masses::Array{T,1},
                              types::Array{Int,1},
                              unitcellmap::AbstractArray;
                              qlinspan::Int = 10,
                              qpoints::Array{S,2} = [ 0 0 0;
                                                      0 0 0],
                              qlabels::Array{Symbol,1} = [:G;
                                                          :G],
                              pdosmesh::Array{Int,1} =[8;8;8],
                              pdosbins::Int = 50) where {T <: Real, S<: Number}

    #Generate reciprocal q-mesh
    rp,tp,qmesh,qnorm = construct_reciprocal_space(cell,
                                                   positions,
                                                   masses,
                                                   qpoints,n=qlinspan)

    #Reciprocal space distance line
    qlabelmap = collect(zip(qlabels,qnorm));
    nqnorm = length(qnorm)-1;                              
    qline = zeros(Float64,nqnorm*qlinspan);
    for i=1:nqnorm
        tmp = range(qnorm[i],stop=qnorm[i+1],length=qlinspan);        
        ii = qlinspan*(i-1)+1:qlinspan*i
        qline[ii] = tmp
    end

    #Get dynamical matrix
    D = getrealspacedynmat(hessian,cell,positions,
                           types,masses);
    Dq = getrecipspacedynmat(D,unitcellmap,qmesh);


    #Eigen Freq. Calc.
    nq, = size(qmesh);
    nd, = size(Dq);
    ν = zeros(Float64,nq,nd);
    for q=1:nq
        ω² = eigvals(Dq[:,:,q]);
        #TODO add check for real valued result
        @assert !(typeof(ω²) <: Complex{T} where {T <: Number}) "ω² at q-$(q) is complex eigenvalue!"
        #TODO Round zero floating point negative values to zero rather than Complex call.
        ω = sort(real.(sqrt.(Complex.(ω²))))
        ν[q,:] = ω .* THz ; 
    end

    #TODO: REWORK AND MOVE
    #Calculate the PDOS using SPGLIB irreducible zone mesh
    #TODO: Assuming Γ-centered meshes, make more general
    pdosqmesh = construct_ir_mesh(cell,
                                  positions,
                                  masses,
                                  pdosmesh)
    
    Dqdos = getrecipspacedynmat(D,unitcellmap,pdosqmesh);
    nq, = size(pdosqmesh);
    nd, = size(Dqdos);
    νdos = zeros(Float64,nq,nd);
    for q=1:nq
        ω² = eigvals(Dqdos[:,:,q]);
        #TODO add check for real valued result
        @assert !(typeof(ω²) <: Complex{T} where {T <: Number}) "ω² at q-$(q) is complex eigenvalue!"
        #TODO Round zero floating point negative values to zero rather than Complex call.
        ω = sort(real.(sqrt.(Complex.(ω²))))
        νdos[q,:] = ω .* THz ; 
    end
 
    #TODO: need to calculate normalization factor
    hfreq,pdos = phonondos(νdos,nbins=pdosbins);

    #Realspace Greens function
    #Local DOS of atom 1
    # eps = 12
    # Ginv = Diagonal(ω².+im*eps) - D ;
    # G = inv(Ginv);
    # ldos = -ω  .* imag(tr(G)) ;
    # open("benchmarks/ldos.txt", "w") do io
    #        writedlm(io,ldos,',')
    #    end;

    return (qlabelmap,qline,ν,hfreq,pdos)
end

function phonons(jobid::String;
                          qlinspan::Int=10,
                          qpoints::Array{S,2}=[0 0 0;
                                               0 0 0],
                          qlabels::Array{Symbol,1} = [:Γ;
                                                      :Γ]) where {S <: Number}
    filename = jobid*".hdf5";
    if isfile(filename) == false
        error("The file $(filename) doesn't exist.")
    end
    
    data = h5open(filename,"r") do file
       read(file,"PotentialCalculation")
    end
    map = read_and_format_unitcellmap_hdf5(filename)

    h = data["Hessian"];
    c = data["cell"]
    p = data["positions"]
    m = data["masses"]
    t = convert(Array{Int},data["types"]);

    #DECISION: create a composite data type for phonon results?
    qlabelmap,qline,ν,hfreq,pdos = phonons(c,h,p,m,t,map,
                              qlinspan=qlinspan,
                              qpoints=qpoints,
                              qlabels=qlabels)

    #Write the phonon dispersion results to the HDF5 file
    writephononhdf5(filename,qpoints,qlabelmap,qline,ν,
                    hfreq,pdos);
    writephonondata(jobid,qpoints,qlabelmap,qline,ν);
    #Write phonon dos, TODO write to HDF5 file.
    writepdos(jobid,hfreq,pdos);
end

@doc """
Get the phonon DOS (histogram) over reciprocal space. Use formula:

``\\(\\rho(E) / V=\\sum_{n} \\sum_{\\mathbf{k}} \\delta\\left(E-E_{n}(\\mathbf{k})\\right)\\)``


""" function phonondos(ν::AbstractArray;
                       nbins::Int=50,
                       normalization::T=1.0) where {T <: Real}
    @info "PDOS calculation not verified."
    nq, = size(ν);
    maxν = maximum(ν);
    minν = 0.00 #minimum(ν);
    diffν = maxν-minν;
    binsize = diffν/nbins;
    histogram = zeros(Float64,nbins+1);
    freqhisto = collect(range(minν,stop=maxν+binsize,length=nbins+1));
    for q=1:nq
        qν = ν[q,:];
        for nu in qν
            ibin = round(Int,(nu*nbins/diffν)) + 1
            histogram[ibin] += 1.0e0
        end
    end
    histogram[1:end] /= normalization;

    return freqhisto,histogram
        
end #phonondos

@doc """
Write the phonon dispersion calculation results to existing HDF5 file.
""" function writephononhdf5(filename::String,
                             qpoints::Array{S,2},
                             qlabelmap::Array{Tuple{Symbol,T},1},
                             qline::Array{Float64,1},
                             ν::Array{Float64,2},
                             hfreq::Array{Float64,1},
                             pdos::Array{Float64,1}) where {S <: Number, T <: Number}
    println("----Writing Phonon Results to HDF5 File-------")
    println(filename)
    h5open(filename,"cw") do file
        #TODO: Need to rework since HDF5.jl throws annoying errors.
        #Delete existing phonon calculation data
        
        try
            o_delete(file,"PhononCalculation");
            @info "Deleting previous PhononCalculation group from HDF5 file."
        catch
            @info "Ingore HDF5 delete error for now!."
        end        
        group = g_create(file,"PhononCalculation");
        group["qlabel"] = [ string(l[1]) for l in qlabelmap];
        attrs(group["qlabel"])["Description"] = "Reciprocal space point labels/symbols"
        group["qlabel-wavenum"] = [float(w[2]) for w in qlabelmap];
        attrs(group["qlabel-wavenum"])["Description"] = "Wavenumber value of reciprocal space labels/symbols"
        attrs(group["qlabel-wavenum"])["Units"] = "Normalized units to Å⁻¹"
        group["qpoints"] = float(qpoints);
        attrs(group["qpoints"])["Description"] = "Reciprocal space points for sampling Brillouin zone" 
        group["wavenumber"] = qline;
        attrs(group["wavenumber"])["Description"] = "Wavenumber normalized length axis.";
        attrs(group["wavenumber"])["Units"] = "Normalized units to Å⁻¹";
        attrs(group["wavenumber"])["Format"] = "N x 1 points"
        nqpoints,nband = size(ν);
        #Loop over each phonon band
        for i=1:nband
            group["band-$i"] = ν[:,i];
            attrs(group["band-$i"])["Description"] = "Phonon dispersion band $(i)";
            attrs(group["band-$i"])["Units"] = "THz"
        end
        #PDOS
        group["pdos-count"] = pdos;
        attrs(group["pdos-count"])["Description"] = "Normalized phonon density of states (PDOS)";
        group["pdos-freq"] = hfreq;
        attrs(group["pdos-freq"])["Description"] = "Bin frequency location of PDOS.";
        attrs(group["pdos-freq"])["Units"] = "THz";
    end
end


@doc """
Write the phonon dispersion calculation results to a data files.
""" function writephonondata(jobid::String,
                             qpoints::Array{S,2},
                             qlabelmap::Array{Tuple{Symbol,T},1},
                             qline::Array{Float64,1},
                             ν::Array{Float64,2}) where {S <: Number, T <: Number}
    bandfile = open(jobid*".bandstructure","w");
    qlabelfile = open(jobid*".bandlabels","w");
    qpointsfile = open(jobid*".qpoints","w");
    
    #Write bandstructure
    npoints,nbands = size(ν);
    write(bandfile,"#q-points $(nbands)-bands \n")
    for n=1:npoints
        qn = qline[n];
        write(bandfile,"$(qn) ");
        for b=1:nbands
            bn=ν[n,b];
            write(bandfile,"$(bn) ");
        end
        write(bandfile,"\n");
    end
    close(bandfile);

    #Write q-labels and locations
    write(qlabelfile,"#q-label norm.-location \n")
    for ql in qlabelmap
        l,p = ql[1],ql[2];
        write(qlabelfile,"$(l) $(p) \n");
    end
    close(qlabelfile)
    
    #Write q-point sampling
    write(qpointsfile,"#q-points sampled in Brillouin zone. \n")
    write(qpointsfile,"#bx by bz \n")
    n, = size(qpoints)
    for qp=1:n
        bx,by,bz = qpoints[qp,1],qpoints[qp,2],qpoints[qp,3];
        write(qpointsfile,"$(bx) $(by) $(bz) \n");
    end
    close(qpointsfile)
end

@doc """
Write the phonon dos to file
""" function writepdos(jobid::String,
                       hfreq::Array,
                       pdos::Array)
    pdosfile = open(jobid*".pdos","w");
    write(pdosfile,"#Bin[THz] Normalized-count \n");
    n = length(hfreq);
    for i=1:n
        b,c = hfreq[i],pdos[i];
        write(pdosfile,"$(b) $(c) \n");
    end
    close(pdosfile)
end #end pdoswrite


    
end #Phonon


       
