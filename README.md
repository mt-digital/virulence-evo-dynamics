# virulence-evo-dnamics

Why do different disease variants with higher or lower virulence (in-host
reproduction rate) evolve? What effect does social structure have on virulence
evolution? This repository implements an agent-based epidemiological 
model with susceptible and infected compartments. Mortality rate and
transmissibility of the model disease are functions of virulence, which evolves
as it is transmitted between agents. The model runs until no more agents are
infected, or up to some user-set maximum time step.


## Installation 

To install and run this project, we use Julia's
[`DrWatson`](https://juliadynamics.github.io/DrWatson.jl/stable/) scientific
project and dependency management toolbox. You will also need a recent R version
installed (the code used here was developed with R v4.2.2. 
There's an open [issue](https://github.com/mt-digital/virulence-evo-dynamics/issues/1) to add R dependency management).

0. Download this code base: `git clone https://github.com/mt-digital/virulence-evo-dynamics`

1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

## Quick start

The following code examples show how to, (1) run a single model parameter
settingin the terminal; (2) run a single model parameter setting in the
terminal, but also use R to visualize the outputs; and (3) run an ensemble
of models with different parameter settings, including multiple computational
trials for each setting.


### Single model run with visualization

Now we need to do the previous codebase `quickactivate` step, but now we want to
`include` code from `scripts/analysis.jl` that contains the `run_series`
function which performs a model run, then plots the susceptible fraction,
infected fraction, 

```julia
using DrWatson; quickactivate("."); include("src/model.jl");
```

```julia
run_series(; whensteps = 10, initial_infected_frac = 0.05, 
             virulence_init = 0.1,  metapop_size = 5000, 
             mutation_rate = 0.8, mutation_variance = 0.2, 
             global_add_rate = 1.5, global_death_rate = 1, 
             maxsteps = 20000, virulence_transmission_coeff=0.02, 
             virulence_mortality_coeff = 0.005, min_homophily = 0.2, 
             maj_homophily = 0.8);
```

This create a plot of the susceptible, infected, and mean virulence for the
total population and broken out for the minority and majority groups. The plot
is saved as `tmp/series.pdf`.


#### Single model run

To see how single model runs work without visualization, see the [`run_series`
source code](scripts/analysis.jl#L9).


### Systematically vary model parameters

For computational analyses we need to run several trials over
systematically-varied model parameters. To do this, Agents.jl provides the
[`ensemblerun!`]() function, which is used in this repository's
`virulence_evo_experiment`, used in the example below. To use this function,
first load the `src/experiment.jl` code.

```julia
using DrWatson; quickactivate("."); include("src/experiment.jl");
```

Then run the experiment using this: 

```julia
resdf = virulence_evo_experiment(10; metapop_size=500, 
                                 virulence_init=collect(0.1:0.2:0.9), 
                                 min_group_frac = 0.2, maxsteps=5000, 
                                 global_add_rate=2, global_death_rate=2, 
                                 min_homophily = 0.8, maj_homophily = 0.2, 
                                 maj_start=true);
```

The output of this can be saved and analysed along with other outputs, say, for
different `virulence_init` values.


## Slurm cluster jobs

[TODO](https://github.com/mt-digital/virulence-evo-dynamics/issues/2)

## Computational analysis of virulence evolution in groups

(WORK IN PROGRESS)

In an R console or R studio, set the project directory as the working directory
(`setwd(".")`) then source the R plotting code: `source("scripts/plot.R")`. The
current draft version of the analysis code is run like so:

```R
draft_total_and_relative_figures()
```

This will put pdf plot outputs in specified folders for the current preliminary 
draft version of this code base. See
[draft_total_and_relative_figures source code](scripts/plot.R#L32).

(TODO: explain how to generate and use output files from computational
experiments)


## Model details

TODO
