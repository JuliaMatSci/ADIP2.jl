

@doc """
Main function to run the potential energy, force, and hessian calculations.
""" function runcalc(jobid::String,
                 potname::Symbol,
                 fullatoms::Array{T,2}) where T <: Real

    #exclude mass column
    atoms = fullatoms[:,2:end];
        
    println("-------------Energy Evaluation------------------")
    energycall = :( $(potname)($(atoms)));
    energy = eval(energycall);
    cohenergy = energy/(length(atoms[:,1]));
    println("------------------------------------------------")
    println("\n")
    println("Total Energy: $(energy) eV")
    println("\n")
    println("Cohesive Energy: $(cohenergy) eV/atom")
    println("\n")
    println("------------------------------------------------")
    println("------------Completed Evaluation----------------")

    #Forces using gradient call using AutoDiff (Forward)
    println("-------------Force Evaluation-------------------")
    expr = :( $(potname) );
    forcecall = a -> -1*ForwardDiff.gradient(eval(expr),a);
    forces = forcecall(atoms)[:,2:end]; #x,y,z
    println("------------Completed Evaluation----------------")

    #TODO: Force constant matrix; gradient call is redundant can be replaced by this call.
    println("-------------Hessian Evaluation-----------------")
    #NOTICE: The atoms attray is .. .. x y z with Natom rows, this means that the
    #calculated Hessian is going to be formated such that it follows the following format
    # dx1-dx1 dx1-dx2 dx1-dxj ... dx1-dy1 dx1-dy2 dx1-dyj ... dx1-dz1 dx1-dz2 dx1-dzj
    # dxi-dxi ................................................................... ...
    # dxi-dyj ................... dyi-dyj ...........................................
    # ..................................................dxi-dzj.............dzi-dzj..
    # It is therefore neccessary to restructure it to be in a more useful format:
    # dxi-dxj dxi-dyj dxi-dzj ...
    # Also need to remove all Zero entries due to the mass and type columns which are
    # meaningless.
    
    #hconfig = ForwardDiff.HessianConfig(eval(potname),atoms,ForwardDiff.Chunk{4}());
    rawhessian = ForwardDiff.hessian(eval(potname),atoms) :: Array{Float64,2}
    hessian = process_raw_hessian(rawhessian,atoms);
    
    println("------------Completed Evaluation----------------")
    
    writehdf5(jobid,potname,energy,fullatoms,forces,hessian)
    
end

using HDF5
@doc """
Write calculation structure, energy, forces, and Hessian to HDF5 file.

DECISION: This maybe better placed in Utilities.jl
""" function writehdf5(jobid::String,
                       potname::Symbol,
                       energy::Float64,
                       atoms::Array{Float64,2},
                       forces::Array{Float64,2},
                       hessian::Array{Float64,2})
    println("-------------Writing HDF5 File-----------------")
    wfile = jobid*".hdf5";
    if isfile(wfile) == false
        error("The file $(wfile) doesn't exist.")
    end
    
    h5open(wfile,"cw") do file
        group = g_create(file,"PotentialCalculation");
        group["potential"] = string(potname);
        group["cell"] = CELL :: Array{Float64,2};
        attrs(group["cell"])["Description"] = "Calculation cell";
        attrs(group["cell"])["Units"] = "Å";
        attrs(group["cell"])["Format"] = "3 X 3 Matrix with coordinate vectors along rows."
        group["energy"] = energy;
        attrs(group["energy"])["Description"] = "Total energy of structure."
        attrs(group["energy"])["Units"] = "eV"
        group["number"] = length(atoms[:,1]);
        attrs(group["number"])["Description"] = "Number of atoms in calculation."
        group["types"] = convert(Array{Int},atoms[:,2]);
        attrs(group["types"])["Description"] = "Indices of atoms corresponding to mass, position, and force arrays."
        group["masses"] = atoms[:,1];
        attrs(group["masses"])["Description"] = "Atomic masses in atomic mass units."
        group["positions"] = atoms[:,3:end];
        attrs(group["positions"])["Description"] = "Atomic coordinates in cartesian representation";
        attrs(group["positions"])["Units"] = "Å"
        attrs(group["positions"])["Format"] = "Rows (x,y,z) X  Columns (Atoms)";
        group["forces"] = forces;
        attrs(group["forces"])["Description"] = "Atomic forces calculated using AutoDiffIP.jl.";
        attrs(group["forces"])["Units"] = "eV/Å";
        attrs(group["forces"])["Format"] = "Rows (fx,fy,fz) X Columns (Atoms)";
        group["Hessian"] = hessian;
        attrs(group["Hessian"])["Description"] = "Atomic Hessian matrix calculated using AutoDiffIP.jl"
        attrs(group["Hessian"])["Units"] = "eV/Å²";
        attrs(group["Hessian"])["Format"] = "3Nₐ X 3Nₐ symmetric matrix";
    end
    
end
