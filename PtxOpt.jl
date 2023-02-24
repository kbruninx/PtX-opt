# Path: PtXOpt.jl/src
module PtxOpt

using PyCall
using CSV
using DataFrames
using Plots
using JSON
using ReusePatterns
using Statistics
using JuMP

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

struct _Parameter_Data
    # Parameter data structure, containing all the parameter data needed for the optimization model
    scenario::Scenario_Data
    components::Dict{String, Union{Component_Data, Electrolyzer_Data}}
    timeseries::Timeseries_Data
end

Base.show(io::IO, p::_Parameter_Data) = print(io, "Parameter Data:
    $(p.scenario),
    $(p.components),
    $(p.timeseries)"
)

# Auxiliary functions
function fetchparameterdata(parameterfile)
    d = JSON.parsefile(parameterfile)
    return _Parameter_Data(
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

function fetchtimeseriesdata(filename::String, scenario_dict::Dict)
    timeseries_df = CSV.read(joinpath(home_dir, filename), DataFrame)
    return Timeseries_Data(
        (gettimeseriesarray(timeseries_df, category, scenario_dict) 
        for category in ["Solar_profile", "Wind_profile", "Dayahead_price"])...
    )
end

function gettimeseriesarray(timeseries_df::DataFrame, category::String, scenario_dict::Dict)
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
#=
    This function takes the hourly values for a whole year and returns a subset of the array
    currently this is just a roughly spaced out subset of days
=#
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
    efficiency::Float64)
    return os_electrolyzer*efficiency*p_electrolyzer
end

function getcomponentcosts(
    capacity::Float64, 
    system_component::Component_Data
)
    return (
        (system_component.crf*system_component.systemprice+system_component.fixopex)*
        timefactor*capacity
    )
end

function getcapitalrecoveryfactor(
    discount_rate::Float64, 
    lifetime::Float64
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



end