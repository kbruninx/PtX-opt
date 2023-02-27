
#= 
The first version of the Julia code for the SET thesis project of Leonard Eble 5357802
Author: Leonard Eble
Date: 2019-01-15
=#

include("PtxOpt.jl")

# Importing the necessary packages
using JuMP
using Gurobi
import .PtxOpt
 

<<<<<<< Updated upstream
function optimizationmodel(parameters::PtxOpt.Parameter_Data)
=======
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

    # Timeseries manipulation
function gettimeseriesarray(
    timeseries_df::DataFrame, 
    category::String, 
    scenario_dict::Dict
)
    hourlyarray = timeseries_df[:, category]
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
        (model[:c_solar], model[:c_wind], model[:c_electrolyzer]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer"))
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
function getelectrolyzerproduction(
    p_electrolyzer::Float64, 
    os_electrolyzer::Bool, 
    efficiency::Float64
)
    return os_electrolyzer*efficiency*p_electrolyzer
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
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
        PtxOpt.lowerbound_electrolyzer(model, electrolyzer, k, t))
=======
        upperbound_electrolyzer(model, electrolyzer, k, t))
>>>>>>> Stashed changes
    ) 

            # Lower bound electrolyzer power
    @constraint(
        model, 
        [k in 1:tc_periods, t in 1:tc_length], 
        (p_electrolyzer[k,t] >= 
        PtxOpt.lowerbound_electrolyzer(model, electrolyzer, k, t))
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
        (sum((PtxOpt.getnetpower(model, electrolyzer, k, t)) for t in 1:tc_length) <= 
        0)
    )

    # Defining the objective function
    @objective(
        model, 
        Min, 
<<<<<<< Updated upstream
        (sum(PtxOpt.getcomponentcosts(parameters, (c_solar, c_wind, c_electrolyzer)) +
        sum(timeseries.price_DAM[k,t]*(p_DAM_buy[k,t]-p_DAM_sell[k,t]) for k in 1:tc_periods, t in 1:tc_length)))
=======
        (getinvestmentcosts(model, parameters) +
        sum(timeseries.price_DAM[k,t]*(p_DAM_buy[k,t]-p_DAM_sell[k,t]) for k in 1:tc_periods, t in 1:tc_length))
>>>>>>> Stashed changes
    )

    # Solving the model
    optimize!(model)

    return model
end

# Printing the results
function showresults(model)
    results = getresultsdict(model)
    # Defining the most commonly used parameternames for readability
    simulation_length = parameters.scenario.simulation_length


<<<<<<< Updated upstream
    # Printing the results
    println("Total costs: ", objective_value(model))
    println("Average cost of Hydrogen: ", objective_value(model)/production_target)
    println("Solar capacity: ", results[:c_solar])
    println("Solar costs: ", getcomponentcosts(results[:c_solar], solar_component))
    println("Wind capacity: ", results[:c_wind])
    println("Wind costs: ", getcomponentcosts(results[:c_wind], wind_component))
    println("Electrolyzer capacity: ", results[:c_electrolyzer])
    println("Electrolyzer costs: ", getcomponentcosts(results[:c_electrolyzer], electrolyzer_component))
    println("CAPEX costs: ", getcomponentcosts(results[:c_solar], solar_component) + getcomponentcosts(results[:c_wind], wind_component) + getcomponentcosts(results[:c_electrolyzer], electrolyzer_component))

    # Preparing the data for plotting
    p_solar = results[:c_solar]*solarcf_subset
    p_wind = results[:c_wind]*windcf_subset
    p_buy = reshape(results[:p_DAM_buy], simulation_length,1)
    p_sell = reshape(results[:p_DAM_sell], simulation_length,1)
    p_el = reshape(results[:p_electrolyzer], simulation_length,1)

    # Printing profit from DAM market
    println("Total profit power market:", sum(price_subset.*(p_sell-p_buy)))
    
    netmarketpowerarray = getnetpower(results[:p_DAM_buy], results[:p_DAM_sell], results[:os_electrolyzer], results[:c_electrolyzer])
    netmarketpower = sum(netmarketpowerarray, dims=2)
    println(netmarketpower)
=======
    capacitycomponentpairs = zip(
        (results[:c_solar], results[:c_wind], results[:c_electrolyzer]),
        (parameters.components[component] for component in ("solar", "wind", "electrolyzer"))
    )

    # Printing the results
    println("Total costs: ", objective_value(model))
    println("Average cost of Hydrogen: ", objective_value(model)/parameters.scenario.production_target)
    
    for (capacity, component) in capacitycomponentpairs
        println("$(name(component)) capacity: ", capacity)
        println("$(name(component)) costs: ", getcomponentcosts(capacity, parameters.scenario, component))
    end

    # Preparing the data for plotting
    p_solar = reshape(results[:c_solar]*parameters.timeseries.solar, simulation_length,1)
    p_wind = reshape(results[:c_wind]*parameters.timeseries.wind, simulation_length,1)
    p_buy = reshape(results[:p_DAM_buy], simulation_length,1)
    p_sell = reshape(results[:p_DAM_sell], simulation_length,1)
    p_el = reshape(results[:p_electrolyzer], simulation_length,1)
    p_DAM = reshape(parameters.timeseries.price_DAM, simulation_length,1)

    # Printing profit from DAM market
    println("Total profit power market:", sum(p_DAM.*(p_sell.-p_buy)))

    legendfontsize = 7
>>>>>>> Stashed changes

    plot(range(1,simulation_length), [p_solar p_wind p_buy p_sell p_el], label=["Solar" "Wind" "Buy" "Sell" "Electrolyzer"], xlabel="Time", ylabel="Power (MW)", title="Power flow", legend=:topleft, legend_font_pointsize=legendfontsize)
    plot!(twinx(), p_DAM, color=:black, linestyle=:dash, label="Price", ylabel="Price (€/MWh)", legend=:topright, show=true, legend_font_pointsize=legendfontsize)
    savefig("powerflow_$(parameters.scenario.name).svg")
end

"""
    main(parameterfilename::String)

"""
function main(parameterfilename::String)
    parameters = PtxOpt.fetchparameterdata(parameterfilename)
    model = optimizationmodel(parameters)
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")
        showresults(model, parameters)
    else
        throw(ErrorException("Optimal solution not found."))
    end
    return model
end
