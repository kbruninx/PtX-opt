
#= 
The first version of the Julia code for the SET thesis project of Leonard Eble 5357802
Author: Leonard Eble
Date: 2023-01-15

The first version consists of:
- RES components (solar and wind)
- Electrolyzer
- DAM market interaction
- Temporal correlation constraints
=#

# Importing the necessary packages
using JuMP
using Gurobi
using PyCall
using CSV
using DataFrames
using Plots
using JSON
using ReusePatterns
using Statistics

# Directory
const home_dir = @__DIR__

# Useful structures
struct Scenario_Data
    name::String
    # Timeframe
    simulation_days::Int64
    simulation_length::Int64 # length in hourly timesteps
    tc_periods::Int64 #number of temporal correlation periods
    tc_length::Int64 #number of hours per temporal correlation period
    # Financial parameters
    discount_rate::Float64
    timefactor::Float64 # Annual to simulation length conversion factor
    # Hydrogen target
    hourly_target::Float64 # kg/h of hydrogen production
    production_target::Float64 # total kg of hydrogen production
end

struct Component_Data
    name:: String
    systemprice::Float64    # €/unit
    fixopex::Float64        # fixed operating expenses as a fraction of the system price
    lifetime::Int64         # lifetime in years
end

# Getters
name(c::Component_Data) = c.name
systemprice(c::Component_Data) = c.systemprice
fixopex(c::Component_Data) = c.fixopex
lifetime(c::Component_Data) = c.lifetime

struct Electrolyzer_Data
    # Electrolyzer data structure, containing additional electrolyzer data needed for the optimization model
    component_data::Component_Data
    # Additional electrolyzer parameters
    efficiency::Float64     # kg/MWh
    P_min::Float64          # MW/MWelz
    P_standby::Float64      # MW/MWelz
end

#Forwards function calls to the component_data field for the Electrolyzer_Data structure
@forward((Electrolyzer_Data, :component_data), Component_Data)

struct Timeseries_Data
    # Timeseries data structure, containing all the timeseries data needed for the optimization model
    solar::Array{Float64,1}
    wind_on::Array{Float64,1}
    wind_off::Array{Float64,1}
    price_DAM::Array{Float64,1}
end

struct Parameter_Data
    # Parameter data structure, containing all the parameter data needed for the optimization model
    scenario::Scenario_Data
    components::Dict{String, Union{Component_Data, Electrolyzer_Data}}
    timeseries::Timeseries_Data
end

#Printing methods
Base.show(io::IO, s::Scenario_Data) = print(io, "Scenario Data:
    Name = $(s.name),
    Simulation days = $(s.simulation_days),
    Simulation length = $(s.simulation_length),
    Temporal Correlation periods = $(s.tc_periods),
    Temporal Correlation length = $(s.tc_length),
    Discount rate = $(s.discount_rate),
    Timefactor = $(s.timefactor),
    Hourly production target = $(s.hourly_target),
    Annual production target = $(s.production_target)" 
)

Base.show(io::IO, c::Component_Data) = print(io, "Component Data:
    Name = $(c.name),
    System price = $(c.systemprice),
    Fixed OPEX = $(c.fixopex),
    Lifetime = $(c.lifetime)"
)

Base.show(io::IO, e::Electrolyzer_Data) = print(io, "Electrolyzer Data:
    $(e.component_data),
    Efficiency = $(e.efficiency),
    Minimal Power = $(e.P_min),
    Standby Power = $(e.P_standby)"
)

Base.show(io::IO, t::Timeseries_Data) = print(io, "Timeseries Data:
    Solar capacity factor mean = $(mean(t.solar)),
    Onshore wind capacity factor mean =  $(mean(t.wind_on)),
    Offshore wind capacity factor mean =  $(mean(t.wind_off)),
    Price DAM mean= $(mean(t.price_DAM))"
)

Base.show(io::IO, p::Parameter_Data) = print(io, "Parameter Data:
    $(p.scenario),
    $(p.components),
    $(p.timeseries)"
)

# Auxiliary functions
    #Fetching data
function fetchparameterdata(parameterfile)
    d = JSON.parsefile(parameterfile)
    return Parameter_Data(
        fetchscenariodata(d["scenario"]),
        fetchcomponentdata(d["components"]),
        fetchtimeseriesdata(d["timeseriesfile"], d["scenario"])
    )
end

function fetchscenariodata(dict)
    return Scenario_Data(
        dict["name"],
        dict["simulation_days"],
        dict["simulation_days"]*24,
        dict["tc_periods"],
        dict["simulation_days"]*24/dict["tc_periods"],
        dict["discount_rate"],
        dict["simulation_days"]*24/(365*24),
        dict["hourly_target"],
        dict["hourly_target"]*dict["simulation_days"]*24
    )
end

function fetchcomponentdata(dict)
    out = Dict()
    for (k, v) in dict
        out[k] = Component_Data(
            k,
            v["systemprice"],
            v["fixopex"],
            v["lifetime"]
        )
    end
    out["electrolyzer"] = Electrolyzer_Data(
        out["electrolyzer"],
        dict["electrolyzer"]["efficiency"],
        dict["electrolyzer"]["P_min"],
        dict["electrolyzer"]["P_standby"] 
    )
    return out
end

# Timeseries manipulation
function fetchtimeseriesdata(
    filename::String, 
    scenario_dict::Dict
)
    timeseries_df = CSV.read(joinpath(home_dir, filename), DataFrame)
    return Timeseries_Data(
        (selectsubset(timeseries_df[:, category], scenario_dict["simulation_days"])
        for category in ["Solar_profile", "Wind_on_profile", "Wind_off_profile", "Dayahead_price"])...
    )
end

function selectsubset(
    hourlyarray::Array{Float64,1}, 
    numberofdays::Int64 
)
    outarray = []
    for i in LinRange(24, length(hourlyarray), numberofdays)
        ni = convert(Int64, round(i))
        append!(outarray, hourlyarray[ni-23:ni])
    end
    return outarray
end

    #Financial calculations

function getcapitalrecoveryfactor(
    discount_rate::Float64, 
    lifetime::Int64
)
    return (discount_rate * (1 + discount_rate)^lifetime / ((1 + discount_rate)^lifetime - 1))
end

function getinvestmentcosts(
    model::JuMP.Model,
    parameters::Parameter_Data,
)
    sum = 0
    capacitycomponentpairs = zip(
        (model[:c_solar], model[:c_wind_on], model[:c_wind_off], model[:c_electrolyzer], model[:c_storage]),
        (parameters.components[component] for component in ("solar", "wind_on", "wind_off", "electrolyzer", "storage"))
    )
    for (capacity, component) in capacitycomponentpairs
        sum += getcomponentcosts(capacity, parameters.scenario, component)
    end
    return sum
end

function getcomponentcosts(
    capacity::Union{VariableRef, Float64}, 
    scenario::Scenario_Data,
    component
)
    crf = getcapitalrecoveryfactor(scenario.discount_rate, lifetime(component))
    return (
        systemprice(component)*capacity*
        scenario.timefactor*
        (crf+fixopex(component))
    )
end

    # Model calculations
function maxcapelectrolyzer(
    parameters::Parameter_Data
)
    min_capacityfactor = 0.05
    return(
        parameters.scenario.hourly_target/parameters.components["electrolyzer"].efficiency/min_capacityfactor
    )
end

function maxcapstorage(
    parameters::Parameter_Data
)
    max_storageoftotal = 0.5
    return(
        max_storageoftotal*parameters.scenario.production_target
    )
end


function getelectrolyzerproduction(
    model::JuMP.Model,
    electrolyzer::Electrolyzer_Data,
    t::Int64
)
    return (
        electrolyzer.efficiency * 
        model[:os_electrolyzer][t] * 
        model[:p_electrolyzer][t]
    )
end

function getnetpower(
    model::JuMP.Model,
    electrolyzer::Electrolyzer_Data,
    t::Int64
)
    return (
        model[:p_DAM_buy][t] .-
        (model[:p_DAM_sell][t] .+
        ((1 .- model[:os_electrolyzer][t]) .*
        (electrolyzer.P_standby .* model[:c_electrolyzer])))
    )
end

function upperbound_electrolyzer(
    model::JuMP.Model, 
    electrolyzer::Electrolyzer_Data, 
    t::Int64
)
    return(
        model[:c_electrolyzer]*model[:os_electrolyzer][t] +
        electrolyzer.P_standby*model[:c_electrolyzer]*(1-model[:os_electrolyzer][t])
    )
end

function lowerbound_electrolyzer(
    model::JuMP.Model, 
    electrolyzer::Electrolyzer_Data, 
    t::Int64
)
    return(
        electrolyzer.P_min*model[:c_electrolyzer]*model[:os_electrolyzer][t] + 
        electrolyzer.P_standby*model[:c_electrolyzer]*(1-model[:os_electrolyzer][t])
    )
end

function getenergyin(
    model::JuMP.Model,
    timeseries::Timeseries_Data,
    t::Int64,
    hourlymatching::Bool
)
    generation = (
        model[:c_solar]*timeseries.solar[t] +
        model[:c_wind_on]*timeseries.wind_on[t] +
        model[:c_wind_off]*timeseries.wind_off[t]
    )
    if hourlymatching
        return generation
    else
        return generation + model[:p_DAM_buy][t]
    end
end

function getmarketvalue(
    model::JuMP.Model,
    timeseries::Timeseries_Data,
    t::Int64,
    hourlymatching::Bool
)
    if hourlymatching
        return -model[:p_DAM_sell][t]*timeseries.price_DAM[t]
    else
        return (model[:p_DAM_buy][t]- model[:p_DAM_sell][t]) * timeseries.price_DAM[t]
    end
end

"""
    

    This function takes the model and returns a dictionary with the results
    The keys are the names of the variables and the values are the values of the variables
"""
function getresultsdict(m::JuMP.Model)
    return Dict(
        k => value.(v) for 
        (k, v) in object_dictionary(m)
    )
end

"""
    solveoptimizationmodel(
        parameters::Parameter_Data
    )

This function takes parameter data and returns a solved model with the results for the given parameters.

# Arguments

- `parameters::Parameter_Data`: Structure containing all the parameters for the optimization model
    
# Returns

- `model::JuMP.Model`: Solved model with the results for the given parameters
"""
function solveoptimizationmodel(
    parameters::Parameter_Data, 
    time_limit::Int64,
    mip_gap::Float64
)
    # Defining the model
    model = Model(Gurobi.Optimizer)

    # Setting the time limit and the mip gap
    set_time_limit_sec(model, time_limit)
    set_optimizer_attribute(model, "MIPgap", mip_gap)

    # Defining the most commonly used parameternames for readability
    simulation_length = parameters.scenario.simulation_length
    tc_length = parameters.scenario.tc_length
    timeseries = parameters.timeseries
    electrolyzer = parameters.components["electrolyzer"]
    hourlymatching = (parameters.scenario.simulation_length == parameters.scenario.tc_periods)
    println("Hourly matching: ", hourlymatching)

    #temporary variables outside JSON files, to be moved later
    S_init = 0.5

    # Defining the variables
        # Design (these caps are currently arbitrary, should be added to json file)
    @variables(model, begin
        0 <= c_solar <= 50 #Solar capacity
        0 <= c_wind_on <= 50 #Onshore capacity
        0 <= c_wind_off <= 50 #Offshore capacity
        0 <= c_electrolyzer <= maxcapelectrolyzer(parameters) #Electrolyzer capacity
        0 <= c_storage <= maxcapstorage(parameters) #Storage capacity
    end)

        # Operation
    @variables(model, begin
        0 <= p_DAM_sell[1:simulation_length] #Power sold to the DAM
        0 <= p_electrolyzer[1:simulation_length] #Electrolyzer power
        os_electrolyzer[1:simulation_length], Bin #Electrolyzer on/standby status
        0 <= h_electrolyzer[1:simulation_length] #Electrolyzer hydrogen production
        0 <= h_storage_soc[1:simulation_length] #Electrolyzer hydrogen storage state of charge
        0 <= h_storage_in[1:simulation_length] #Electrolyzer hydrogen storage input
        0 <= h_storage_out[1:simulation_length] #Electrolyzer hydrogen storage output
        0 <= h_demand[1:simulation_length] #Hydrogen output
    end)

    #Currently, it is really slow for hourly simulations, this is a temporary fix, where buying is excluded for hourly
    if !(hourlymatching)

            # Power bought from the DAM
        @variable(
            model,
            0 <= p_DAM_buy[1:simulation_length]
        )
        
            # Prevents buying and selling at the same time, speeding up model
        @constraint(
            model,
            [t in 1:simulation_length],
            p_DAM_buy[t]*p_DAM_sell[t] == 0
        )
    
            # Upper bound power bought from the DAM, which should never exceed the max rated power of the electrolyzer
        @constraint(
            model, 
            [t in 1:simulation_length], 
            p_DAM_buy[t] <= c_electrolyzer
        )

            #Temporal correlation constraint
        @constraint(
            model, 
            [k in (1:tc_length:simulation_length)], 
            (sum((getnetpower(model, electrolyzer,t)) for t in (k:(k+tc_length-1))) <= 
            0)
        )
    end

    # Defining the constraints
        # Number of upper bounds on variables that dramatically increase the model speed
            # Upper bound electrolyzer capacity, which should never exceed the max rated power of the RES
    @constraint(
        model, 
        (c_electrolyzer <= (c_wind_on+c_wind_off+c_solar))
    )
    
            # Upper bound power sold to the DAM, which could never exceed the max rated power of the RES
    @constraint(
        model, 
        [t in 1:simulation_length],
        p_DAM_sell[t] <= (c_wind_on+c_wind_off+c_solar)
    ) 

            # Upper bound electrolyzer hydrogen storage input, which should never exceed the max rated power of the electrolyzer
    @constraint(
        model,
        [t in 1:simulation_length],
        h_storage_in[t] <= parameters.scenario.hourly_target*5 #temporary fix, should be replaced by c_compressor or something of the sort
    )

    @constraint(
        model,
        [t in 1:simulation_length],
        h_storage_out[t] <= parameters.scenario.hourly_target*5
    )

        # Physical constraints
            # Energy balance
    @constraint(
        model, 
        [t in 1:simulation_length], 
        (p_electrolyzer[t] + p_DAM_sell[t] == 
        getenergyin(model, timeseries, t, hourlymatching))
    )

            # Upper bound electrolyzer power
    @constraint(
        model, 
        [t in 1:simulation_length], 
        (p_electrolyzer[t] <= 
        upperbound_electrolyzer(model, electrolyzer, t))
    ) 

            # Lower bound electrolyzer power
    @constraint(
        model, 
        [t in 1:simulation_length], 
        (p_electrolyzer[t] >= 
        lowerbound_electrolyzer(model, electrolyzer, t))
    )

        # Electrolyzer hydrogen production
    @constraint(
        model,
        [t in 1:simulation_length],
        (h_electrolyzer[t] == getelectrolyzerproduction(model, electrolyzer, t))
    )
    # Hydrogen storage capacity
    @constraint(
        model,
        [t in 1:simulation_length],
        (h_storage_soc[t] <= c_storage)
    )
    # Hydrogen output
    @constraint(
        model,
        [t in 1:simulation_length],
        (h_demand[t] == h_electrolyzer[t] + h_storage_out[t] - h_storage_in[t])
    )
        # Hydrogen storage initialization
    @constraint(
        model,
        (h_storage_soc[1] == S_init*c_storage + h_storage_in[1]*0.95 - h_storage_out[1])
    )
        # Electrolyzer hydrogen storage state of charge
    @constraint(
        model,
        [t in 2:simulation_length],
        (h_storage_soc[t] == h_storage_soc[t-1] + h_storage_in[t]*0.95 - h_storage_out[t])
    )
        # Minimal hourly production
    @constraint(
        model,
        [t in 1:simulation_length],
        h_demand[t] >= 0.8*parameters.scenario.hourly_target
    )
        # Production target
    @constraint(
        model, 
        (sum(h_demand[t] for t in 1:simulation_length) >=
        (parameters.scenario.production_target+S_init*c_storage))
    ) 

    # Defining the objective function
    @objective(
        model, 
        Min, 
        (getinvestmentcosts(model, parameters) +
        sum(getmarketvalue(model, timeseries, t, hourlymatching) for t in 1:simulation_length))
    )

    # Solving the model
    optimize!(model)

    return model
end

"""
    getresults(
        model::JuMP.Model, 
        parameters::Parameter_Data,
        verbose::Bool
    )

This function returns the results of the optimization model.

# Arguments

- `model::JuMP.Model`: Solved model with the results for the given parameters
- `parameters::Parameter_Data`: Structure containing all the parameters for the optimization model

# Returns

- `var_data::Dict{Symbol,Array{Float64,1}}`: Dictionary containing the results of the optimization model
"""
function getresults(
    model::JuMP.Model, 
    parameters::Parameter_Data,
    verbose::Bool
)
    var_data = getresultsdict(model)
    # Defining the most commonly used parameternames for readability
    simulation_length = parameters.scenario.simulation_length

    capacitycomponentpairs = zip(
        (var_data[:c_solar], var_data[:c_wind_on], var_data[:c_wind_off], var_data[:c_electrolyzer], var_data[:c_storage]),
        (parameters.components[component] for component in ("solar", "wind_on", "wind_off", "electrolyzer", "storage"))
    )

    # Preparing the data for plotting
    p_solar = parameters.timeseries.solar*var_data[:c_solar]
    p_wind_on = parameters.timeseries.wind_on*var_data[:c_wind_on]
    p_wind_off = parameters.timeseries.wind_off*var_data[:c_wind_off]
    p_DAM = parameters.timeseries.price_DAM
    p_DAM_buy = 0
    p_DAM_sell = var_data[:p_DAM_sell]

    outcome_data = DataFrame(name = String[], capacity= Any[], cost = Any[])

    push!(outcome_data, ("Total costs", missing, objective_value(model)))
    push!(outcome_data, ("Average cost of Hydrogen", missing, objective_value(model)/parameters.scenario.production_target))
    push!(outcome_data, ("Electrolyzer capacity factor", mean(var_data[:p_electrolyzer]./var_data[:c_electrolyzer]), missing))
    for (capacity, component) in capacitycomponentpairs
        push!(outcome_data, (name(component), capacity, getcomponentcosts(capacity, parameters.scenario, component)))
    end
    push!(outcome_data, ("Total profit power market", missing, sum(p_DAM.*(p_DAM_sell.-p_DAM_buy))))

    # Printing the results if verbose
    if verbose
        print(outcome_data)
    end

    legendfontsize = 3

    # Collecting all the plotted data
    powerlabels = ["Solar", "Wind (onshore)", "Wind (offshore)", "Sell", "Electrolyzer"]
    powerflow = [p_solar, p_wind_on, p_wind_off, var_data[:p_DAM_sell], var_data[:p_electrolyzer]]
    hydrogenflow = [var_data[:h_electrolyzer], var_data[:h_storage_in], var_data[:h_storage_out], var_data[:h_demand]]

    if haskey(var_data, :p_DAM_buy)
        println("Including buy from DAM market")
        p_DAM_buy = var_data[:p_DAM_buy]
        insert!(powerlabels, 4, "Buy")
        insert!(powerflow, 4, p_DAM_buy)
    end

    #transforms labels back into matrix, which is needed for plotting
    powerlabels = permutedims(powerlabels)

    # Plotting the results and adding them to the plot dict
    powerplot = plot(
        range(1,simulation_length), powerflow, 
        label=powerlabels, 
        xlabel="Time", ylabel="Power (MW)", title="Power flow", 
        legend=:topleft, legend_font_pointsize=legendfontsize)
    # Plots the price of the DAM market on a second y-axis
    powerplot = plot!(
        twinx(), p_DAM, color=:black, linestyle=:dash, 
        label="Price", ylabel="Price (€/MWh)", legend=:topright, 
        legend_font_pointsize=legendfontsize)
    
    
    hydrogenplot = plot(
        range(1, simulation_length), hydrogenflow, 
        label=["Electrolyzer" "Storage in" "Storage out" "Demand"], 
        xlabel="Time", ylabel="Hydrogen (kg)", title="Hydrogen flow", legend=:topleft, 
        legend_font_pointsize=legendfontsize)
    # Plots the storage SOC on a second y-axis
    hydrogenplot = plot!(
        twinx(), var_data[:h_storage_soc]./var_data[:c_storage], color=:black, linestyle=:dash, 
        label="Storage SOC", ylabel="Storage SOC (%)", legend=:topright, 
        legend_font_pointsize=legendfontsize)
    plots = Dict("powerflow" => powerplot, "hydrogenflow" => hydrogenplot)

    #Plot subsets if simulation is long
    if simulation_length >= 24*14
        subsetsize = 24*7
        subsetpowerflow = map(x->x[1:subsetsize], powerflow)
        subsethydrogenflow = map(x->x[1:subsetsize], hydrogenflow)

        sub_powerplot = plot(
            range(1, subsetsize), subsetpowerflow, 
            label=powerlabels, 
            xlabel="Time", ylabel="Power (MW)", title="Power flow", 
            legend=:topleft, legend_font_pointsize=legendfontsize)
        sub_powerplot = plot!(
            twinx(), p_DAM[1:subsetsize], color=:black, 
            linestyle=:dash, label="Price", ylabel="Price (€/MWh)", 
            legend=:topright, legend_font_pointsize=legendfontsize)

        sub_hydrogenplot = plot(
            range(1, subsetsize), subsethydrogenflow, 
            label=["Electrolyzer" "Storage in" "Storage out" "Demand"], 
            xlabel="Time", ylabel="Hydrogen (kg)", 
            title="Hydrogen flow", legend=:topleft, 
            legend_font_pointsize=legendfontsize)
        sub_hydrogenplot = plot!(
            twinx(), var_data[:h_storage_soc][1:subsetsize]./var_data[:c_storage], 
            color=:black, linestyle=:dash, label="Storage SOC", 
            ylabel="Storage SOC (%)", legend=:topright, 
            legend_font_pointsize=legendfontsize)
        subplots = Dict("powerflow(subplot)" => sub_powerplot, "hydrogenflow(subplot)" => sub_hydrogenplot)
        plots = merge(plots, subplots)
    end

    return Dict("data" => outcome_data, "figures" => plots)
end

function saveresults(
    results::Dict{String, Any},
    parameters::Parameter_Data,
    savepath::String,
    savecsv::Bool,
    savefigs::Bool
)
    if savecsv
        CSV.write("$(savepath)/variabledata_$(parameters.scenario.name).csv", results["data"])
    end
    if savefigs
        for (k,v) in results["figures"]
            savefig(v, "$(savepath)/$(k)_$(parameters.scenario.name).svg")
        end
    end
end

"""
    main(parameterfilename::String)

This function is the main function of the optimization model. 
    It takes a parameterfilename and solves the optimization model with the given parameters.

# Arguments

- `parameterfilename::String`: Name of the parameterfile, which should be a .json file

# Returns

- `model::JuMP.Model`: Solved model with the results for the given parameters

"""
function main(
    parameterfilename::String, 
    verbose::Bool=true,
    savepath::String="$(home_dir)/Figures",
    savefigs::Bool=true,
    savecsv::Bool=true
)
    parameters = fetchparameterdata(parameterfilename)
    if verbose
        println("Parameters: ", parameters)
    end
    model = solveoptimizationmodel(parameters, 240, 0.05)
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")
        results = getresults(model, parameters, verbose)
        saveresults(results, parameters, savepath, savefigs, savecsv)
    else
        throw(ErrorException("Optimal solution not found."))
    end
    return results
end