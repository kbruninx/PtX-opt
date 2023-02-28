
#=
Testing file for V1.jl
=#
# Directory
const home_dir = @__DIR__

using JuMP
using Gurobi
using PyCall
using CSV
using DataFrames
using Plots
using JSON
using ReusePatterns
using Statistics

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
    solar::Array{Float64,2}
    wind::Array{Float64,2}
    price_DAM::Array{Float64,2}
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

# TESTING
using Test

testscenario1 = Scenario_Data("test1", 200, 8760, 8760, 1, 0.05, 8760/8760, 1, 8760)

testtimeseriesrand = Timeseries_Data(rand(4,5), rand(4,5), rand(4,5))
testtimeseries1 = Timeseries_Data(
    [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5],
    [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5],
    [1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5; 1 2 3 4 5]
)

testcomponent1 = Component_Data("test1", 1000, 100, 20)
testcomponent2  = Component_Data("test2", 2000, 200, 30)
testcomponent3 = Electrolyzer_Data(Component_Data("test2", 2000, 200, 30), 0.8, 0.1, 0.2)
testcomponents = Dict("test1" => testcomponent1, "test2" => testcomponent2, "test3" => testcomponent3)

testparameters1 = Parameter_Data(testscenario1, testcomponents, testtimeseries1)

testcapacity1 = 100
testcapacity2 = 200
testcapacity3 = 300

capacities = Dict("solar" => testcapacity1, "wind" => testcapacity2, "electrolyzer" => testcapacity3)


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
        (model[:c_solar], model[:c_wind], model[:c_electrolyzer]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer"))
    )
    for (capacity, component) in capacitycomponentpairs
        sum += getcomponentcosts(capacity, component, parameters.scenario)
    end
    return sum
end

function getcomponentcosts(
    capacity::Union{VariableRef, Float64, Int64}, 
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


# function confirmedcomponentcosts(
#     capacity::Float64,
#     component::Component_Data,
#     scenario::Scenario_Data
# )
#     crf = getcapitalrecoveryfactor(scenario.discount_rate, scenario.life_time)
#     cost = capacity * scenario.timefactor * (crf*systemprice(component) + fixopex(component))
#     return cost
# end

@testset "getcomponentcosts" begin
    @test getcomponentcosts(testcapacity1, testscenario1, testcomponent1) == 1000
    @test getcomponentcosts(testcapacity2, testscenario1, testcomponent2) == 2000
    @test getcomponentcosts(testcapacity3, testscenario1, testcomponent3) == 2000
end