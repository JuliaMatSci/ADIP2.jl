# Author: Stefan Bringuier
# Email: stefanbringuier@gmail.com
# see LICENSE file for usage details.

# Info: This is the main file to launch the program ADIP2.jl
# The file executes the following:
# 1. Read the input Julia script given as a command-line argument
# 2. Parse the potential info and construct potential call.
# 3. Execute routines for eneergy and ForwardDiff.jl.
# 4. Perform phonon analysis if requested.


#TODO: convert this into a parser function for text input file.
# for example the goal would be to use a LAMMPS input file to
# specify configuration and potential/parameter file
include(ARGS[1])


#TODO: make automatic
#SET EXE DIRECTORY
BASE="/home/stefanb/GoogleDrive/COMPUTATION/ADIP2.jl/src/"
#cd(BASE)
push!(LOAD_PATH,BASE)


#External Packages
using LinearAlgebra
   #AutoDiff Packages
using ForwardDiff
#using Zygote
#using ReverseDiff


#Internal Modules
  using Potentials
  using GenerateParameters
  using Utilities
  using PrepareCalc
  using DynamicalMatrix
  using Phonon



printbanner()


#### CHECK VARIABLE EXISTANCE ####
if ! @isdefined jobid
    jobid="ADIP_Calculation";
end




if ! @isdefined basis
    error("A variable named `` basis `` is REQUIRED!")
end


#### GET POTENTIAL DETAILS ####
if ! @isdefined potinfo
    error("A variable named `` potinfo `` is REQUIRED!")
else
    potname = potinfo[:potname];
    potentialdriver = potentialdatabase()[potname];
    include(BASE*potentialdriver)
end
### DEFINE IAP REQUIRED CONSTANTS ####
global const PARAMETERS = genallparams(potinfo)
###############################

#SKIP TO PHONON if flagged
pexist = @isdefined phonononly
if !pexist || phonononly == false

    ### GENERATE COMPUTATION ATOMIC STRUCTURE ###
    if ! @isdefined(cell)
        error("A variable named `` cell `` is REQUIRED!")
    elseif ! @isdefined(replicate)
        error("A variable named `` replicate `` is REQUIRED!")
    else


        #TODO: WOULD LIKE TO PUT ELSEWHERE ###
        #Get primitive cell along with unit cell mapping
        #and  reconstituted ``` atoms``` array.
        bpos = basis[:,3:end];
        bmass = basis[:,1];
        btypes = convert(Array{Int},basis[:,2])    
        if ! @isdefined(spglib)
            unitcellmap, extendedcell, atoms = generate_atomic_structure(cell,
                                                             bpos,
                                                             bmass,
                                                             btypes,
                                                         replicate=replicate);
        elseif @isdefined(spglib) && (spglib == false)
            unitcellmap, extendedcell, atoms = generate_atomic_structure(cell,
                                                             bpos,
                                                             bmass,
                                                             btypes,
                                                             replicate=replicate,
                                                             spglib=false);
        end
       write_unitcellmap_hdf5(jobid,unitcellmap);
    end

    ### DEFINE CONSTANT CELL WHICH IS REQUIRED TO USE AutoDiff ###
    global const CELL = extendedcell
    #######################################

    #### RUN PROGRAM ####
    if ! @isdefined atoms
        error("A variable named `` atoms `` is REQUIRED!")
    else    
        #MAIN FILE with functions using ForwardDiff.jl
        include(BASE*"RunCalc.jl")
        runcalc(jobid,potname,atoms)
    end
    #### END PROGRAM ####
end

#### RUN ANALYSIS ####
#TODO - rework using metaprogramming construct
if @isdefined phonon
    if phonon == true
        println("--------Calculating Phonon Dispersion-----------")
        if (@isdefined qpoints) && (@isdefined qlinspan)
            @info "Assuming that ```qlabels``` variable exist"
            phonons(jobid,qlinspan=qlinspan,
                             qpoints=qpoints, qlabels=qlabels);
        elseif (@isdefined qpoints) && !(@isdefined qlinspan)
            @info "Assuming that ```qlabels``` variable exist"
            phonons(jobid,qpoints=qpoints,qlabels=qlabels);
        elseif !(@isdefined qpoints) && (@isdefined qlinspan)
            phonons(jobid,qlinspan=qlinspan);
        else
            phonondispersion(jobid);
        end
        println("------------Calculation Complete----------------")
    end
end

#### END ANALYSIS #### 
