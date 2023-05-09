# PtX-opt
PtX component optimization - MSc project Leonard

This project aims to calculate optimal hydrogen production costs, given temporal constraints on both electricity sourcing and hydrogen delivery. 
It optimizes for rated capacities (design) and operations using an MILP model based on perfect foresight. 

List of files:

* OptimizationModel.jl

The main module which runs the optimization given a number of parameters. The function mainscript() takes as a main input a .JSON filepath, in the format of the parameterdata_Template.JSON file, which contains all the parameters needed to run the optimization. 
It returns the results, with hourly data on hydrogen flows and power flows, and an outcome dictionary, with design and cost outcomes. This is either returned as a single dictionary, or in the case of a scenario analysis in the form of a dataframe. 
It is currently recommended to calculate any desired results and save them to a variable, which can then be fed to functions such as makeplots() or makeheatmap() to generate visual representations. 

* parameterdata_Template.JSON

Contains all the parameter information. Parameter data can be divided in three categories: scenario data, component data, and timeseries data. 
For data references, see the Maths.TeX file.
To run a number of different scenario's, the inputs for the scenario data can be given in the form of arrays (e.g., no_etc_periods=[12,365,8760]). In this case, each permutation of each of the input values is calculated.

* ES90_–.csv

Example data drawn from the ENTSOE transparency platform for 01.2020-01.2023. All data is for the bidding zone ES. The capacity factors are calculated by normalizing the 
actual generation using the maximum value for the year. This could eventually be adjusted, considering this does not adequately account for new generation being installed throughout the year, causing the values earlier in the year to be underestimated.
The data has been transformed to a set of representative days with decision weightings and orderings representing each non-represented day as a linear-representation of represented days.

* csvtools.jl

A collection of functions that have been used to edit data taken from the ENTSOE transparency platform, which could be useful for future csv maniupulation.

* runtests.jl

Unittesting file intended for OptimizationModel.jl. Currently not functional, due to problems with submodules in Julia.

* Project.toml | Manifest.toml 

Files used by julia to ensure the correct packages are installed.