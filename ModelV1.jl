#= 
The first version of the Julia code for the SET thesis project of Leonard Eble 5357802
Author: Leonard Eble
Date: 2019-01-15
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

Base.show(io::IO, s::Scenario_Data) = print(io, "Scenario Data:
    Simulation days = $(s.simulation_days),
    Simulation length = $(s.simulation_length),
    Temporal Correlation periods = $(s.tc_periods),
    Temporal Correlation length = $(s.tc_length),
    Discount rate = $(s.discount_rate),
    Timefactor = $(s.timefactor),
    Hourly production target = $(s.hourly_target),
    Annual production target = $(s.production_target)" 
)

struct Component_Data
    name:: String
    systemprice::Float64    # €/unit
    fixopex::Float64        # fixed operating expenses as a fraction of the system price
    lifetime::Int64         # lifetime in years
end

systemprice(c::Component_Data) = c.systemprice
fixopex(c::Component_Data) = c.fixopex
lifetime(c::Component_Data) = c.lifetime

Base.show(io::IO, c::Component_Data) = print(io, "Component Data:
    Name = $(c.name),
    System price = $(c.systemprice),
    Fixed OPEX = $(c.fixopex),
    Lifetime = $(c.lifetime)"
)

struct Electrolyzer_Data
    # Electrolyzer data structure, containing additional electrolyzer data needed for the optimization model
    component_data::Component_Data
    # Additional electrolyzer parameters
    efficiency::Float64     # kg/MWh
    P_min::Float64          # MW/MWelz
    P_standby::Float64      # MW/MWelz
end

@forward((Electrolyzer_Data, :component_data), Component_Data)

Base.show(io::IO, e::Electrolyzer_Data) = print(io, "Electrolyzer Data:
    $(e.component_data),
    Efficiency = $(e.efficiency),
    Minimal Power = $(e.P_min),
    Standby Power = $(e.P_standby)"
)

struct Timeseries_Data
    # Timeseries data structure, containing all the timeseries data needed for the optimization model
    solar::Array{Float64,2}
    wind::Array{Float64,2}
    price_DAM::Array{Float64,2}
end

Base.show(io::IO, t::Timeseries_Data) = print(io, "Timeseries Data:
    Solar capacity factor mean = $(mean(t.solar)),
    Wind capacity factor mean =  $(mean(t.wind)),
    Price DAM mean= $(mean(t.price_DAM))"
)

struct Parameter_Data
    # Parameter data structure, containing all the parameter data needed for the optimization model
    scenario::Scenario_Data
    components::Dict{String, Union{Component_Data, Electrolyzer_Data}}
    timeseries::Timeseries_Data
end

Base.show(io::IO, p::Parameter_Data) = print(io, "Parameter Data:
    $(p.scenario),
    $(p.components),
    $(p.timeseries)"
)

# Auxiliary functions
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

function fetchtimeseriesdata(
    filename::String, 
    scenario_dict::Dict
)
    timeseries_df = CSV.read(joinpath(home_dir, filename), DataFrame)
    return Timeseries_Data(
        (gettimeseriesarray(timeseries_df, category, scenario_dict) 
        for category in ["Solar_profile", "Wind_profile", "Dayahead_price"])...
    )
end

function gettimeseriesarray(
    timeseries_df::DataFrame, 
    category::String, 
    scenario_dict::Dict
)
    hourlyarray = timeseries_df[:, category]
    #println(
    #    "Working on $category: $(typeof(hourlyarray))
    #    The problem is at $(findall(x->typeof(x)==Missing, hourlyarray))")
    subset = selectsubset(hourlyarray, scenario_dict["simulation_days"])
    return reshape(
        subset, 
        scenario_dict["tc_periods"], 
        convert(Int64, scenario_dict["simulation_days"]*24/scenario_dict["tc_periods"])
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

function getelectrolyzerproduction(
    p_electrolyzer::Float64, 
    os_electrolyzer::Bool, 
    efficiency::Float64
)
    return os_electrolyzer*efficiency*p_electrolyzer
end

function getcomponentcosts(
    model::JuMP.Model,
    parameters::Parameter_Data,
)
    sum = 0
    capacitycomponentpairs = zip(
        (model[:c_solar], model[:c_wind], model[:c_electrolyzer]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer"))
    )
    for (capacity, component) in capacitycomponentpairs
        crf = getcapitalrecoveryfactor(
            parameters.scenario.discount_rate, 
            lifetime(component)
        )
        sum += (
            crf*systemprice(component)+
            fixopex(component))*parameters.scenario.timefactor*capacity
    end
    return sum
end

function getcapitalrecoveryfactor(
    discount_rate::Float64, 
    lifetime::Int64
)
    return (discount_rate * (1 + discount_rate)^lifetime / ((1 + discount_rate)^lifetime - 1))
end

function getnetpower(
    model::JuMP.Model,
    electrolyzer::Electrolyzer_Data,
    k::Int64, t::Int64
)
    return (
        model[:p_DAM_buy][k,t] .-
        (model[:p_DAM_sell][k,t] .+
        ((1 .- model[:os_electrolyzer][k,t]) .*
        (electrolyzer.P_standby .* model[:c_electrolyzer])))
    )
end

function upperbound_electrolyzer(
    model::JuMP.Model, 
    electrolyzer::Electrolyzer_Data, 
    k::Int64, t::Int64
)
    return(
        model[:c_electrolyzer]*model[:os_electrolyzer][k,t] +
        electrolyzer.P_standby*model[:c_electrolyzer]*(1-model[:os_electrolyzer][k,t])
    )
end

function lowerbound_electrolyzer(
    model::JuMP.Model, 
    electrolyzer::Electrolyzer_Data, 
    k::Int64, t::Int64
)
    return(
        electrolyzer.P_min*model[:c_electrolyzer]*model[:os_electrolyzer][k,t] + 
        electrolyzer.P_standby*model[:c_electrolyzer]*(1-model[:os_electrolyzer][k,t])
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



function optimizationmodel(parameters::Parameter_Data)
    # Defining the model
    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, 60)

    # Defining the most commonly used parameternames for readability
    tc_periods = parameters.scenario.tc_periods
    tc_length = parameters.scenario.tc_length
    timeseries = parameters.timeseries
    electrolyzer = parameters.components["electrolyzer"]


    # Defining the variables
        # Design (POTENTIALLY ADD UPPER BOUNDS HERE)
    @variables(model, begin
        0 <= c_solar #Solar capacity
        0 <= c_wind #Wind capacity
        0 <= c_electrolyzer #Electrolyzer capacity
    end)

        # Operation
    @variables(model, begin
        0 <= p_DAM_buy[1:tc_periods, 1:tc_length] #Power bought from the DAM
        0 <= p_DAM_sell[1:tc_periods, 1:tc_length] #Power sold to the DAM
        0 <= p_electrolyzer[1:tc_periods, 1:tc_length] #Electrolyzer power
        os_electrolyzer[1:tc_periods, 1:tc_length], Bin #Electrolyzer on/standby status
    end)

    # Defining the constraints
        # Number of upper bounds on variables that dramatically increase the model speed
            # Upper bound electrolyzer capacity, which should never exceed the max rated power of the RES
    @constraint(
        model, 
        (c_electrolyzer <= (c_wind+c_solar))
    )

            # Upper bound power bought from the DAM, which should never exceed the max rated power of the electrolyzer
    @constraint(
        model, 
        [k in 1:tc_periods, t in 1:tc_length], 
        p_DAM_buy[k,t] <= c_electrolyzer
    )
    
            # Upper bound power sold to the DAM, which could never exceed the max rated power of the RES
    @constraint(
        model, 
        [k in 1:tc_periods, t in 1:tc_length],
        p_DAM_sell[k,t] <= (c_wind+c_solar)
    ) 

        # Physical constraints
            # Energy balance
    @constraint(
        model, 
        [k in 1:tc_periods, t in 1:tc_length], 
        (p_electrolyzer[k,t] + p_DAM_sell[k,t] == 
        c_solar*timeseries.solar[k,t] + c_wind*timeseries.wind[k,t]+ p_DAM_buy[k,t])
    )

            # Upper bound electrolyzer power
    @constraint(
        model, 
        [k in 1:tc_periods, t in 1:tc_length], 
        (p_electrolyzer[k,t] <= 
        lowerbound_electrolyzer(model, electrolyzer, k, t))
    ) 

            # Lower bound electrolyzer power
    @constraint(
        model, 
        [k in 1:tc_periods, t in 1:tc_length], 
        (p_electrolyzer[k,t] >= 
        lowerbound_electrolyzer(model, electrolyzer, k, t))
    )

        # Production target
    @constraint(
        model, 
        (sum(os_electrolyzer[k,t]*electrolyzer.efficiency*p_electrolyzer[k,t] for k in 1:tc_periods, t in 1:tc_length) >=
        parameters.scenario.production_target)
    )

        # Temporal correlation
    @constraint(
        model, 
        [k in 1:tc_periods], 
        (sum((getnetpower(model, electrolyzer, k, t)) for t in 1:tc_length) <= 
        0)
    )

    # Defining the objective function
    @objective(
        model, 
        Min, 
        (sum(getcomponentcosts(model, parameters) +
        sum(timeseries.price_DAM[k,t]*(p_DAM_buy[k,t]-p_DAM_sell[k,t]) for k in 1:tc_periods, t in 1:tc_length)))
    )

    # Solving the model
    optimize!(model)

    return model
end

# Printing the results
function showresults(
    model::JuMP.Model, 
    parameters::Parameter_Data
)
    results = getresultsdict(model)

    capacitycomponentpairs = zip(
        (model[:c_solar], model[:c_wind], model[:c_electrolyzer]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer"))
    )

    # Printing the results
    println("Total costs: ", objective_value(model))
    println("Average cost of Hydrogen: ", objective_value(model)/parameters.scenario.production_target)
    
    for (capacity, component) in capacitycomponentpairs
        println("$(component.name) capacity: ", capacity)
        println("$(component.name) costs: NONE")
    end

    # Preparing the data for plotting
    p_solar = results[:c_solar]*parameters.timeseries.solar
    p_wind = results[:c_wind]*parameters.timeseries.wind
    p_buy = reshape(results[:p_DAM_buy], parameters.scenario.simulation_length,1)
    p_sell = reshape(results[:p_DAM_sell], parameters.scenario.simulation_length,1)
    p_el = reshape(results[:p_electrolyzer], parameters.scenario.simulation_length,1)

    # Printing profit from DAM market
    println("Total profit power market:", sum(parameters.timeseries.price_DAM.*(p_sell-p_buy)))

    plot(range(1,simulation_length), [p_solar p_wind p_buy p_sell p_el], label=["Solar" "Wind" "Buy" "Sell" "Electrolyzer"], xlabel="Time", ylabel="Power (MW)", title="Power flow", legend=:topleft, legend_font_pointsize=legendfontsize)
    plot!(twinx(), price_subset, color=:black, linestyle=:dash, label="Price", ylabel="Price (€/MWh)", legend=:topright, show=true, legend_font_pointsize=legendfontsize)
    savefig("powerflow.svg")
end

"""
    main(parameterfilename::String)

"""
function main(parameterfilename::String)
    parameters = fetchparameterdata(parameterfilename)
    model = optimizationmodel(parameters)
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")
        showresults(model)
    else
        throw(ErrorException("Optimal solution not found."))
    end
    return model
end