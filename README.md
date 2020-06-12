# MATLAB Traveling Salesman Problem (TSP) Genetic Algorithm Toolbox

## Purpose
This toolbox contains MATLAB functions to solve the Traveling Salesman Problem (TSP), Multiple Traveling Salesman Problem (M-TSP) and other variations using a custom Genetic Algorithm (GA)

[![View Traveling Salesman Problem (TSP) Genetic Algorithm Toolbox on the MATLAB File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/75525-traveling-salesman-problem-tsp-genetic-algorithm-toolbox)

## TSP/M-TSP and Variations
Functions in this toolbox are listed below with short descriptions

* `tsp_ga` - Solves the classic Traveling Salesman Problem (TSP)
* `tspo_ga` - Solves an open variation of the TSP with no start/end constraints
* `tspof_ga` - Solves an open variation of the TSP with fixed start & end constraints
* `tspofs_ga` - Solves an open variation of the TSP with only a fixed start constraint
* `tsp_ga_max` - Solves a variation of the TSP that *maximizes* the total distance
* `tsp_nn` - Nearest Neighbor (NN) solution to the TSP
* `tspo_nn` - Nearest Neighbor (NN) solution to the open variation of the TSP (no start/end constraints)
* `tsp_fn` - Farthest Neighbor (FN) solution to the TSP
* `tsp_rs` - Random Search (RS) solution to the TSP (provides interesting performance comparisons)
* `tsp_ga_co` - Solves the TSP with a variation of the GA that uses a cross-over operator
* `tsp_ga_hybrid` - Solves the TSP with a variation of the GA that uses a hybrid set of operators
* `tsp_ga_clusters` - Variation of the TSP to find the shortest tour through clusters of points
* `tsp_ga_turbo` - Solves the TSP with a variation of the GA that increases the mutation rate on the best route
* `tspo_ga_turbo` - Solves the Open TSP with a variation of the GA that increases the mutation rate on the best route
* `seed_tsp` - Script to compare `tsp_ga` performance with an approach that seeds the algorithm with the NN solution
* `seed_tsp_turbo` - Script to compare `tsp_ga_turbo` performance with an approach that seeds the algorithm with the NN solution
* `mtsp_ga` - Solves the classic Multiple Traveling Salesmen Problem (M-TSP)
* `mtsp_rs` - Random Search (RS) solution to the M-TSP (provides interesting performance comparisons)
* `mtspf_ga` - Solves a variation of the M-TSP where all salesmen start and end at the first city
* `mtspo_ga` - Solves an open variation of the M-TSP with no start/end constraints
* `mtspof_ga` - Solves an open variation of the M-TSP where all salesmen start at the first city and end at the last
* `mtspofs_ga` - Solves an open variation of the M-TSP where all salesmen start at the first city
* `mtspf_ga_bases` - Solves a variation of the M-TSP where all salesmen are assigned starting cities (bases)
* `mtspof_ga_bases` - Solves an open variation of the M-TSP where all salesmen are assigned separate starting and ending cities (bases)
* `mtspofs_ga_bases` - Solves an open variation of the M-TSP where all salesmen are assigned separate starting cities (bases)
* `mtspofs_ga_depots` - Solves an open variation of the M-TSP where all salesmen are assigned separate starting cities and must travel to a specified end city (depot) that is not assigned in advance
* `mtsp_ga_max` - Solves a variation of the M-TSP that *maximizes* the total distance
* `mtsp_ga_minsum` - Same as `mtsp_ga` (just named such that objective function is explicit)
* `mtsp_ga_minmax` - Solves the classic M-TSP except the objective function minimizes the maximum tour (which tends to make tour lengths more equitable)
* `mtsp_ga_combo` - Solves the classic M-TSP except the objective function is a combination of *minmax* and *minsum*
* `mtspf_ga_minmax` - Solves the `mtspf_ga` variation, but with the *minmax* objective function
* `mtsp_ga_turbo` - Solves the M-TSP with a variation of the GA that increases the mutation rate on the best route
* `mtspv_ga` - Solves a variation of `mtsp_ga` that also optimizes the number of salesmen
* `mtspvf_ga` - Solves a variation of `mtspf_ga` that also optimizes the number of salesmen
* `mtspvo_ga` - Solves a variation of `mtspo_ga` that also optimizes the number of salesmen
* `mtspvof_ga` - Solves a variation of `mtspof_ga` that also optimizes the number of salesmen
* `mtspvofs_ga` - Solves a variation of `mtspofs_ga` that also optimizes the number of salesmen
* `figstatus` - Helper function to display a progress status bar along the bottom of the figure with the option to stop early
* `get_config` - Helper function to support a variety of input options (structure or name/value pairs)

## Notes about Custom GA
The custom Genetic Algorithm used by most of the functions in this toolbox does not use crossover and mutation operators in the traditional sense, because the crossover operator tends to be a highly destructive operator and rarely improves the best solution. Instead, three different mutation operators (*flip*, *swap*, and *slide*) are used -- see `TSPGA-Mutation-Descriptions.pdf` for more details.

## Features
The solvers in this toolbox have several (hopefully useful) features:

1. The inputs can either be a structure with zero or more expected fields, or a set of name/value pairs
2. The output structure from running any of the GA solvers can be passed back as an input into the solver to continue where it left off
3. If the plot showing the current best solution is displayed, a status bar is shown across the bottom of the figure with an option to exit early
4. If the waitbar is used, it also utilizes an option that lets the user exit early if desired
5. The output structures from the solvers contain fields ("plotPoints", "plotResult", "plotMatrix", and "plotHistory") that hold anonymous functions which can be used as a simple way to view the results (see examples)

## Examples
Below are examples demonstrating some of the toolbox features

Solve the classic TSP with a randomly generated set of city locations and then inject the result structure back into the solver as an input to continue where it left off

	cities = 10*rand(100,2);
    tsp = tsp_ga('xy',cities,'numIter',1000);
	tsp = tsp_ga(tsp);

Solve the classic M-TSP without any plots/status and then use the anonymous function fields to show the results

	cities = 10*rand(30,2);
    mtsp = mtsp_ga('xy',cities,'showProg',false,'showResult',false);
	figure; mtsp.plotPoints(mtsp); title('Points')
	figure; mtsp.plotMatrix(mtsp); title('Sorted Distance Matrix')
	figure; hold('on'); mtsp.plotResult(mtsp); title('Solution')
	figure; mtsp.plotHistory(mtsp); title('Distance History')

