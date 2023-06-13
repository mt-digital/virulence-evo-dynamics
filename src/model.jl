##
# Agent-based SIR model of infectious disease spread where the pathogen's
# virulence (represented here by recovery rate) evolves to an optimal value
# for the given group structure.

# I drew on Simon Frost's Agents.jl SIR model as a helpful example 
# (https://github.com/epirecipes/sir-julia/blob/master/script/abm/abm.jl)
# for the SIR portion of this model. 

# I combined an SIR model with code I wrote for understanding the spread of
# adaptations in metapopulations, currently at 
# https://github.com/eehh-stanford/SustainableCBA/blob/main/src/model.jl. 
# But here we have disease transmission instead of social learning as 
# in SustainbleCBA.

# Author: Matthew A. Turner <maturner01@gmail.com>
# Date: May 25, 2023
#

using Agents
using Distributions
using Random
using StatsBase


"""
Classic SIR infectious disease model states–we leave Dead agents in simulation for performance reasons.
"""
@enum SIR_Status Susceptible Infected #Dead


"""
Pathogen agents only hold their virulence; we assume a virulence-transmissibility
tradeoff where increased transmissibility comes from infected agents being 
infected longer (lower recovery rate), at the risk of increased mortality
rate for the pathogen–the mortality rate is calculated using the heuristic that
virulence + mortality_rate = 1, so if the recovery rate is 0.6 the mortality
rate is 0.4, i.e., there is a 60% chance an infected Person recovers on a 
given time step, and if they don't recover, there is a 40% chance the 
infected Person dies. 
"""
struct Pathogen
    # recovery_rate::Float64
    virulence::Float64
end


@enum Group Minority Majority


mutable struct Person <: AbstractAgent
    
    id::Int

    status::SIR_Status

    group::Group

    # Each agent is infected by a unique pathogen with a certain virulence.
    pathogen::Pathogen

    homophily::Float64
end


"""
Agent steps will use people; pathogens evolve on transmission,
within this agent_step!, for simplicity.
"""
function agent_step!(focal_agent::Person, model::ABM)

    virulence = copy(focal_agent.pathogen.virulence)

    # Possibly get infected if not infected...
    if focal_agent.status == Susceptible
        interact!(focal_agent, model)

    # ...or possibly die if infected
    elseif (focal_agent.status == Infected) &&
           (rand() < virulence_mortality_rate(
                        virulence; 
                        virulence_mortality_coeff = model.virulence_mortality_coeff
                    )
           )
        
        remove_agent!(focal_agent, model)
    end
end


function model_step!(model)

    # Possibly add agents through birth/migration.
    add_count = rand(model.add_count_dist)

    if add_count > 0
        for _ in 1:add_count
            if rand() < model.min_group_frac
                add_agent!(model, Susceptible, Minority, Pathogen(NaN), 
                           model.min_homophily)
            else
                add_agent!(model, Susceptible, Majority, Pathogen(NaN), 
                           model.maj_homophily)
            end
        end
    end

    # Possibly remove random agents.
    death_count = rand(model.death_count_dist)

    if death_count > 0
        for _ in 1:death_count
            random_agent_dieoff!(model)
        end
    end
end


function random_agent_dieoff!(model)

    agent_to_remove = sample(collect(allagents(model)))
    remove_agent!(agent_to_remove, model)

    # agent_to_remove.status = Dead
end


"""
Virulence tradeoff: more virulent means lower recovery rate, higher death rate.
"""
function virulence_mortality_rate(virulence::Float64; 
                                  virulence_mortality_coeff = 0.05)

    # For now use linear virulence-mortality relationship.
    return virulence_mortality_coeff * virulence
end


function sample_group(focal_agent, model)

    weights = zeros(2)

    # XXX a waste to calculate this every time.
    agent_group_weight = (1 + focal_agent.homophily) / 2.0

    group_idx = Int(focal_agent.group) + 1

    # weights[focal_agent.group] = agent_group_weight
    # weights[1:end .!= focal_agent.group] .= 1 - agent_group_weight
    weights[group_idx] = agent_group_weight
    weights[1:end .!= group_idx] .= 1 - agent_group_weight
    
    return sample([Minority, Majority], Weights(weights)) 
end


"""
Select interaction partner, interact, possibly get infected, and if infection
happens, let the pathogen evolve.
"""
function interact!(focal_agent, model)

    # partner = sample(allagents(model))
    group = sample_group(focal_agent, model)
    partner = select_partner(focal_agent, model, group)
    while partner == focal_agent
        partner = sample(allagents(model))
    end

    # Possibly get infected.
    if (partner.status == Infected) && 
       (rand() ≤ transmissibility(partner.pathogen.virulence; 
                                  a = model.virulence_transmission_coeff,
                                  b = model.virulence_transmission_denom_summand))

        focal_agent.status = Infected 

        # Without mutation, the infection has the same recovery rate as partner's.
        virulence = copy(partner.pathogen.virulence)

        # The pathogen evolves.
        if rand() < model.mutation_rate
            virulence += rand(model.mutation_dist)

            if virulence < 0.0
                virulence = 0.05
            end
        end

        # Initialize new Pathogen object for the transmitted pathogen.
        transmitted_pathogen = Pathogen(virulence)
        focal_agent.pathogen = transmitted_pathogen

        model.total_infected += 1
        if focal_agent.group == Minority
            model.total_minority_infected += 1
        else
            model.total_majority_infected += 1
        end
    end
end


function virulence_evo_model(; metapop_size = 1000, virulence_init = 0.3,
                               initial_infected_frac = 0.10, mutation_rate = 0.05,
                               mutation_variance = 0.05, 
                               global_add_rate = 0.0, global_death_rate = 0.0,
                               virulence_mortality_coeff = 0.1, 
                               virulence_transmission_coeff = 1.0, 
                               virulence_transmission_denom_summand = 5.0,
                               min_group_frac = 0.2,
                               min_start = true,
                               maj_start = false,
                               min_homophily = 0.0,
                               maj_homophily = 0.0,
                               rep_idx = nothing)

    # Mutations are drawn from normal distros with zero mean and given variance.
    mutation_dist = Normal(0.0, mutation_variance)

    # Death count for each time step drawn from Poisson with λ = death_rate.
    add_count_dist = Poisson(global_add_rate)

    # Death count for each time step drawn from Poisson with λ = death_rate.
    death_count_dist = Poisson(global_death_rate)

    # Track total number of infections over time.
    total_infected::Int = 0
    total_minority_infected::Int = 0
    total_majority_infected::Int = 0

    properties = @dict(metapop_size, min_group_frac, min_start, maj_start,
                       min_homophily, maj_homophily,
                       virulence_init, initial_infected_frac, 
                       virulence_mortality_coeff,
                       virulence_transmission_coeff,
                       virulence_transmission_denom_summand,
                       mutation_rate, mutation_dist, global_death_rate, 
                       global_add_rate, death_count_dist, add_count_dist,
                       total_infected, total_minority_infected,
                       total_majority_infected, rep_idx)

    # model = UnremovableABM(Person; properties)
    model = ABM(Person; properties)

    initialize_metapopulation!(model)
    
    return model
end


"""
Given a minority fraction; population size; group-level homophily values; which
group starts with the infection, and the ABM where the agents will live, create
and add the appropriate number of agents from each group to the population.
"""
function initialize_metapopulation!(model::ABM)
        
    virulence_init = model.virulence_init
    metapop_size = model.metapop_size
    initial_infected_frac = model.initial_infected_frac

    virulence_init_distro = Normal(virulence_init, 0.1)

    gmin_count = Int(model.min_group_frac * metapop_size)
    gmaj_count = metapop_size - gmin_count

    if model.min_start
        gmin_init_infected_count = ceil(initial_infected_frac * gmin_count)
    else
        gmin_init_infected_count = 0
    end
    for min_agent_idx in 1:gmin_count
        initialize_add_agent!(min_agent_idx, model, gmin_init_infected_count,
                             Minority)
    end

    if model.maj_start
        gmaj_init_infected_count = ceil(initial_infected_frac * gmaj_count)
    else
        gmaj_init_infected_count = 0
    end
    for maj_agent_idx in 1:gmaj_count
        initialize_add_agent!(maj_agent_idx, model, 
                              gmaj_init_infected_count,
                              Majority, gmin_count)
    end

end



function initialize_add_agent!(in_group_idx, model, initial_infected_count, 
                               group, agent_model_index_offset = 0)

    if in_group_idx ≤ initial_infected_count

        status = Infected
        # virulence_init = rand(virulence_init_distro)
        virulence_init = model.virulence_init

        if virulence_init < 0.0
            virulence_init = 0.0
        elseif virulence_init > 1.0
            virulence_init = 1.0
        end

        pathogen = Pathogen(virulence_init)
    else

        status = Susceptible
        pathogen = Pathogen(NaN)
    end

    model_agent_idx = in_group_idx + agent_model_index_offset

    homophily = group == Minority ? model.min_homophily : model.maj_homophily

    add_agent!(
        Person(model_agent_idx, status, group, pathogen, homophily), model
    )
end


# Classic saturating transmissibility-virulence relationship.
function transmissibility(virulence; a = 1.0, b = 10.0)
    return (a * virulence) / (b + virulence)
end


function select_partner(focal_agent, model, group)

    ## Begin payoff-biased social learning from teacher within selected group.
    prospective_partners = 
        filter(agent -> (agent.group == group) && 
               (agent != focal_agent),# &&
               # (agent.status != Dead), 
               collect(allagents(model)))

    # partner_weights = 
    #     map(agent -> model.trait_fitness_dict[agent.curr_trait], 
    #                           prospective_partners)

    # # Renormalize weights.
    # denom = Float64(sum(partner_weights))
    # partner_weights ./= denom

    # Select partner.
    return sample(prospective_partners)  #, Weights(partner_weights))
end
