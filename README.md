# PtX-opt
PtX component optimization - MSc project Leonard

This project aims to calculate optimal hydrogen production costs, given temporal constraints on both electricity sourcing and hydrogen delivery. 
It optimizes for rated capacities (design) and operations using an MILP model based on perfect foresight. 

List of files:

* OptimizationModel.jl

The main module which runs the optimization given a number of parameters. The function main() takes as a main input the parameterdata.JSON file, 
which contains all the parameters needed to run the optimization. It returns the results, which it also saves to csv files and figures, given that 
the function is run with the savefigs and savecsv arguments as true. These files are saved in the Results folder.

* parameterdata.JSON

Contains all the parameter information. Parameter data can be divided in three categories: scenario data, component data, and timeseries data. 
For data references, see the Maths.TeX file.

* Timeseriesdata 2021.csv

Example data drawn from the ENTSOE transparency platform. All data is for the bidding zone NL. The capacity factors are calculated by normalizing the 
actual generation using the maximum value for the year. This should eventually be adjusted, considering this does not adequately account for new generation
being installed throughout the year, causing the values earlier in the year to be underestimated.

* csvtools.jl

A collection of functions that have been used to edit data taken from the ENTSOE transparency platform, which could be useful for future csv maniupulation.

* runtests.jl

Unittesting file intended for OptimizationModel.jl. Currently not functional, due to problems with submodules in Julia.

* Project.toml | Manifest.toml 

Files used by julia to ensure the correct packages are installed.