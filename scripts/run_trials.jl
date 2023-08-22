using Distributed
using Dates
using UUIDs: uuid4

using DrWatson
# quickactivate("..")
quickactivate(".")

# Set up DrWatson to include vectors in autogen save names.
da = DrWatson.default_allowed(Any)
DrWatson.default_allowed(c) = (da..., Vector)

using ArgParse
using CSV
using Comonicon
using JLD2

include("../src/experiment.jl")


s = ArgParseSettings()


function parse_cli()

    @add_arg_table s begin
        "datadirname"
            help = "Subdirectory of data/, where to save output data"
            arg_type = String
            required = true

        "--nreplicates"
            help = "Number of trial simulations to run for this experiment"
            arg_type = Int
            default = 100

        "--metapop_size", "-N"
            help = "Population size of all groups combined, N"
            default = 2000
            arg_type = Int

        "--min_group_frac", "-m"
            help = "Fraction of population that is in minority"
            arg_type = Float64
            default = 0.05

        "--min_start"
            help = "Whether the minority group has infected agents"
            arg_type = Bool
            default = true

        "--maj_start"
            help = "Whether the majority group has infected agents"
            arg_type = Bool
            default = true

        "--virulence_init"
            help = "Initial virulence shared by all initially-infected agents"
            default = collect(0.1:0.2:0.9)
            # required = true

        "--initial_infected_frac"
            help = "Initial fraction in each group who are infected"
            default = 0.05

        "--virulence_mortality_coeff"
            help = "Linear coefficient specifying how mortality rate scales with virulence"
            default = 0.005

        "--mutation_rate"
            help = "Virulence mutation rate"
            default = 0.8

        "--mutation_variance"
            help = "Virulence mutation variance"
            default = 0.2

        "--virulence_transmission_coeff"
            help = "Multiplicative coefficient in numerator of transmissibility as a function of virulence"
            default = 0.02

        "--min_homophily"
            help = "Minority group homophily level"
            default = 0.0

        "--maj_homophily"
            help = "Majority group homophily level"
            default = 0.0

        "--global_add_rate"
            help = "Birth rate/rate of migration in to metapopulation"
            default = 1.5

        "--global_death_rate"
            help = "Death rate/rate of migration out of metapopulation"
            default = 1.0

    end

    return parse_args(s)
end



function run_trials(nreplicates = 20; 
                    outputfilename = "trials_output.jld2", 
                    experiment_kwargs...)

    tic = now()

    println("Starting trials at $(replace(string(tic), "T" => " "))")

    # XXX Awkward stuff due to mixing around positional argument as either
    # nagents or nreplicates.
    kwargs_dict = Dict(experiment_kwargs)
    # nagents = pop!(kwargs_dict, :nagents)
    # kwargs_dict[:nreplicates] = nreplicates

    result_df = virulence_evo_experiment(nreplicates; kwargs_dict...)

    CSV.write(outputfilename, result_df)

    trialstime = Dates.toms(now() - tic) / (60.0 * 1000.0)

    println("Ran expected payoffs trials in $trialstime minutes")
end


function main()
    parsed_args = parse_cli()

    # Create job id for unique filename.
    parsed_args["jobid"] = string(uuid4())
    println(parsed_args)

    datadirname = pop!(parsed_args, "datadirname")
    nameargs = copy(parsed_args)

    outputfilename = joinpath("data", datadirname, savename(nameargs, "csv"))

    nreplicates = pop!(parsed_args, "nreplicates")

    # Don't need to pass this job ID to experiment.
    pop!(parsed_args, "jobid")

    # Need keys to be symbols for passing to run_trials function.
    pa_symbkeys = Dict(Symbol(key) => value for (key, value) in parsed_args)

    run_trials(nreplicates; outputfilename, pa_symbkeys...)
end


main()
