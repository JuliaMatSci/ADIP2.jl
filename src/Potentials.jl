#push!(LOAD_PATH,".")
__precompile__()
#List of all potential files with call functions

module Potentials
export potentialdatabase

@doc """
Return a single dictionary with the potential keys and associated string values.
The string value is the name of the potential routine, e.g., CallTersoff.jl. Simply add 
any new potentials """ function potentialdatabase() 

    Dict(:tersoff => "CallTersoff.jl",
         :lennardjones => "CallLennardJones.jl",
         :eamsutton => "CallEAM.jl");
    end
end
