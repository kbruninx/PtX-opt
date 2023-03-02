
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
    wind::Array{Float64,1}
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
    Wind capacity factor mean =  $(mean(t.wind)),
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
        for category in ["Solar_profile", "Wind_profile", "Dayahead_price"])...
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

function getTCsets(
    scenario::Scenario_Data
)   
    simlength = scenario.simulation_length
    tclength = scenario.tc_length
    return (
        [(i:(i+tclength-1)) for i in (1:tclength:simlength)]
    )
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
        (model[:c_solar], model[:c_wind], model[:c_electrolyzer], model[:c_storage]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer", "storage"))
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
function solveoptimizationmodel(parameters::Parameter_Data)
    # Defining the model
    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, 120)
    set_optimizer_attribute(model, "MIPgap", 0.05)

    # Defining the most commonly used parameternames for readability
    simulation_length = parameters.scenario.simulation_length
    tc_length = parameters.scenario.tc_length
    timeseries = parameters.timeseries
    electrolyzer = parameters.components["electrolyzer"]

    #temporary variables outside json
    S_init = 0.5

    # Defining the variables
        # Design (POTENTIALLY ADD UPPER BOUNDS HERE)
    @variables(model, begin
        0 <= c_solar <= 50 #Solar capacity
        0 <= c_wind <= 50#Wind capacity
        0 <= c_electrolyzer <= maxcapelectrolyzer(parameters) #Electrolyzer capacity
        0 <= c_storage <= maxcapstorage(parameters) #Storage capacity
    end)

        # Operation
    @variables(model, begin
        0 <= p_DAM_buy[1:simulation_length] #Power bought from the DAM
        0 <= p_DAM_sell[1:simulation_length] #Power sold to the DAM
        0 <= p_electrolyzer[1:simulation_length] #Electrolyzer power
        os_electrolyzer[1:simulation_length], Bin #Electrolyzer on/standby status
        0 <= h_electrolyzer[1:simulation_length] #Electrolyzer hydrogen production
        0 <= h_storage_soc[1:simulation_length] #Electrolyzer hydrogen storage state of charge
        0 <= h_storage_in[1:simulation_length] #Electrolyzer hydrogen storage input
        0 <= h_storage_out[1:simulation_length] #Electrolyzer hydrogen storage output
        0 <= h_demand[1:simulation_length] #Hydrogen output
    end)

    # Defining the constraints
        # Number of upper bounds on variables that dramatically increase the model speed
            # Upper bound electrolyzer capacity, which should never exceed the max rated power of the RES
    @constraint(
        model,
        [t in 1:simulation_length],
        p_DAM_buy[t]*p_DAM_sell[t] == 0
    )

    @constraint(
        model, 
        (c_electrolyzer <= (c_wind+c_solar))
    )

            # Upper bound power bought from the DAM, which should never exceed the max rated power of the electrolyzer
    @constraint(
        model, 
        [t in 1:simulation_length], 
        p_DAM_buy[t] <= c_electrolyzer
    )
    
            # Upper bound power sold to the DAM, which could never exceed the max rated power of the RES
    @constraint(
        model, 
        [t in 1:simulation_length],
        p_DAM_sell[t] <= (c_wind+c_solar)
    ) 

            # Upper bound electrolyzer hydrogen storage input, which should never exceed the max rated power of the electrolyzer
    @constraint(
        model,
        [t in 1:simulation_length],
        h_storage_in[t] <= 10000
    )

    @constraint(
        model,
        [t in 1:simulation_length],
        h_storage_out[t] <= 10000
    )

        # Physical constraints
            # Energy balance
    @constraint(
        model, 
        [t in 1:simulation_length], 
        (p_electrolyzer[t] + p_DAM_sell[t] == 
        c_solar*timeseries.solar[t] + c_wind*timeseries.wind[t]+ p_DAM_buy[t])
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
        h_demand[t] >= parameters.scenario.hourly_target
    )
        # Production target
    @constraint(
        model, 
        (sum(h_demand[t] for t in 1:simulation_length) >=
        (parameters.scenario.production_target+S_init*c_storage))
    )

    @constraint(
        model, 
        [k in (1:tc_length:simulation_length)], 
        (sum((getnetpower(model, electrolyzer,t)) for t in (k:(k+tc_length-1))) <= 
        0)
    )

    # Defining the objective function
    @objective(
        model, 
        Min, 
        (getinvestmentcosts(model, parameters) +
        sum(timeseries.price_DAM[t]*(p_DAM_buy[t]-p_DAM_sell[t]) for t in 1:simulation_length))
    )

    # Solving the model
    optimize!(model)

    return model
end

"""
    showresults(
        model::JuMP.Model, 
        parameters::Parameter_Data
    )

This function takes a solved model and the parameters and prints a set of results for the optimization model.
    It also saves the powerflow results to a csv file.

# Arguments

- `model::JuMP.Model`: Solved model with the results for the given parameters
- `parameters::Parameter_Data`: Structure containing all the parameters for the optimization model

# Returns

- `nothing`
"""
function showresults(
    model::JuMP.Model, 
    parameters::Parameter_Data
)
    results = getresultsdict(model)
    # Defining the most commonly used parameternames for readability
    simulation_length = parameters.scenario.simulation_length

    capacitycomponentpairs = zip(
        (results[:c_solar], results[:c_wind], results[:c_electrolyzer], results[:c_storage]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer", "storage"))
    )

    # Printing the results
    println("Total costs: ", objective_value(model))
    println("Average cost of Hydrogen: ", objective_value(model)/parameters.scenario.production_target)
    
    componentcosts = Dict()
    for (capacity, component) in capacitycomponentpairs
        componentcosts[name(component)] = getcomponentcosts(capacity, parameters.scenario, component)
        println("$(name(component)) capacity: ", capacity)
        println("$(name(component)) costs: ", componentcosts[name(component)])
    end

    # Preparing the data for plotting
    p_solar = parameters.timeseries.solar*results[:c_solar]
    p_wind = parameters.timeseries.wind*results[:c_wind]
    println("LCOE is: ", (sum(p_solar.+p_wind)/(componentcosts["solar"]+componentcosts["wind"])))
    powerflow = [p_solar, p_wind, results[:p_DAM_buy], results[:p_DAM_sell], results[:p_electrolyzer]]
    hydrogenflow = [results[:h_electrolyzer], results[:h_storage_in], results[:h_storage_out], results[:h_demand]]
    dfh = DataFrame(hydrogenflow, :auto)
    println(first(dfh, 5))
    p_DAM = parameters.timeseries.price_DAM
    # Printing profit from DAM market
    println("Total profit power market:", sum(p_DAM.*(results[:p_DAM_sell].-results[:p_DAM_buy])))

    legendfontsize = 7

    plot(range(1,simulation_length), powerflow, label=["Solar" "Wind" "Buy" "Sell" "Electrolyzer"], xlabel="Time", ylabel="Power (MW)", title="Power flow", legend=:topleft, legend_font_pointsize=legendfontsize)
    plot!(twinx(), p_DAM, color=:black, linestyle=:dash, label="Price", ylabel="Price (€/MWh)", legend=:topright, show=true, legend_font_pointsize=legendfontsize)
    savefig("powerflow_$(parameters.scenario.name).svg")
    plot(range(1, simulation_length), hydrogenflow, label=["Electrolyzer" "Storage in" "Storage out" "Demand"], xlabel="Time", ylabel="Hydrogen (kg)", title="Hydrogen flow", legend=:topleft, legend_font_pointsize=legendfontsize)
    savefig("hydrogenflow_$(parameters.scenario.name).svg")
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
function main(parameterfilename::String)
    parameters = fetchparameterdata(parameterfilename)
    model = solveoptimizationmodel(parameters)
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")
        showresults(model, parameters)
    else
        throw(ErrorException("Optimal solution not found."))
    end
    return model
end


