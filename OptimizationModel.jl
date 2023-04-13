
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
using JSON
using ReusePatterns
using Statistics
using CairoMakie
CairoMakie.activate!(type = "svg")

# Directory
const home_dir = @__DIR__

# Useful structures
struct Scenario_Data
    name::String
    # Timeframe
    simulation_days::Int64
    simulation_length::Int64 # length in hourly timesteps
    # Financial parameters
    discount_rate::Float64
    lambda_TSO::Float64 # TSO penalty factor
    timefactor::Float64 # Annual to simulation length conversion factor
    # Hydrogen target
    hourly_target::Float64 # kg/h of hydrogen production
    production_target::Float64 # kg of hydrogen production
    etc_length::Int64 #number of hours per electricity temporal correlation period
    htc_length::Int64 #number of hours per hyrdogen temporal correlation period
    flexrange::Float64 #flexibility range in % of the hourly target
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

struct Storage_Data
    # Storage data structure, containing additional storage data needed for the optimization model
    component_data::Component_Data
    # Additional storage parameters
    efficiency::Float64     # kgout/kgin
    soc_init::Float64       # %
end

#Forwards function calls to the component_data field for the Storage_Data structure
@forward((Storage_Data, :component_data), Component_Data)

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
    components::Dict{String, Union{Component_Data, Electrolyzer_Data, Storage_Data}}
    timeseries::Timeseries_Data
end

#Printing methods
Base.show(io::IO, s::Scenario_Data) = print(io, "Scenario Data:
    Name = $(s.name),
    Simulation days = $(s.simulation_days),
    Simulation length = $(s.simulation_length),
    Discount rate = $(s.discount_rate),
    TSO tariff = $(s.lambda_TSO),
    Timefactor = $(s.timefactor),
    Hourly production target = $(s.hourly_target),
    Total production target = $(s.production_target),
    Electrical Temporal Correlation length = $(s.etc_length),
    Hydrogen Temporal Correlation length = $(s.htc_length),
    Flexibility range = $(s.flexrange)"
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

Base.show(io::IO, s::Storage_Data) = print(io, "Storage Data:
    $(s.component_data),
    Efficiency = $(s.efficiency),
    Initial State of Charge = $(s.soc_init)"
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
    scenarios = combinescenarios(d["scenarios"])
    components = componentsfromdict(d["components"])
    timeseries = timeseriesfromdict(d["timeseriesfile"], d["scenarios"])

    # If there is only one scenario, return a Parameter_Data
    if scenarios isa Scenario_Data
        return Parameter_Data(
                scenarios,
                components,
                timeseries
        )
    # If there are multiple scenarios, return a DataFrame with a column of Parameter_Data
    elseif scenarios isa DataFrame
        parameterdata = select(scenarios, Not("Scenario"))
        parameterdata[:, "Parameters"] = [Parameter_Data(
            scenario,
            components,
            timeseries
            ) for scenario in scenarios[:, "Scenario"]]
        return parameterdata
    else
        error("Scenario data is neither a Scenario_Data nor a DataFrame")
    end
end

function permutatedict(
    dict::Dict{String, Any}
)
    #Checking for variables of interest (arrays)
    varsofinterest = [key for (key, value) in dict if value isa Array]
    #If there are none, return the scenario
    if varsofinterest == []
        return scenariofromdict(dict), varsofinterest
    end
    #Take out the name of scenario and permutate the rest
    name = pop!(dict, "name")
    permutations = [Dict{String,Any}(zip(keys(dict), p)) for p in Iterators.product(values(dict)...)]
    #Put the name back in
    for p in permutations
        p["name"] = name
    end

    return permutations, varsofinterest
end

function combinescenarios(
    scenariodict::Dict{String, Any}
    )
    #Get the scenarios and the variables of interest
    scenarios, varsofinterest = permutatedict(scenariodict)

    #If there are no variables of interest, return a Scenario_Data
    if varsofinterest == []
        return scenarios
    end

    #Initialize the DataFrame
    df = DataFrame([columnname => [] for columnname in [varsofinterest..., "Scenario"]])
    #Fill the DataFrame
    for p in scenarios
        push!(df, ([p[var] for var in varsofinterest]..., scenariofromdict(p)))
    end
    return df
end

function scenariofromdict(dict)
    return Scenario_Data(
        dict["name"],
        dict["simulation_days"],
        dict["simulation_days"]*24,
        dict["discount_rate"],
        dict["lambda_TSO"],
        dict["simulation_days"]*24/(365*24),
        dict["hourly_target"],
        dict["hourly_target"]*dict["simulation_days"]*24,
        dict["etc_length"],
        dict["htc_length"],
        dict["flexrange"]
    )
end

function componentsfromdict(dict)
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
    out["storage"] = Storage_Data(
        out["storage"],
        dict["storage"]["efficiency"],
        dict["storage"]["soc_init"]
    )
    return out
end

# Timeseries manipulation
function timeseriesfromdict(
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
    min_capacityfactor = 0.05 #REVIEW this value
    return(
        parameters.scenario.hourly_target/parameters.components["electrolyzer"].efficiency/min_capacityfactor
    )
end

function maxcapstorage(
    parameters::Parameter_Data
)
    max_storageoftotal = 0.5 #REVIEW this value
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
    t::Int64
)
    generation = (
        model[:c_solar]*timeseries.solar[t] +
        model[:c_wind_on]*timeseries.wind_on[t] +
        model[:c_wind_off]*timeseries.wind_off[t]
    )
    return generation + model[:p_DAM_buy][t]
end

function getmarketvalue(
    model::JuMP.Model,
    timeseries::Timeseries_Data,
    lambda_TSO::Float64,
    t::Int64
)
    return model[:p_DAM_buy][t] * (timeseries.price_DAM[t] + lambda_TSO) - model[:p_DAM_sell][t] * timeseries.price_DAM[t]
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
    set_optimizer_attribute(model, "DualReductions", 0)

    # Defining the most commonly used parameternames for readability
    simulation_length = parameters.scenario.simulation_length
    hourly_target = parameters.scenario.hourly_target
    etc_length = parameters.scenario.etc_length
    htc_length = parameters.scenario.htc_length
    flexrange = parameters.scenario.flexrange
    timeseries = parameters.timeseries
    lambda_TSO = parameters.scenario.lambda_TSO
    electrolyzer = parameters.components["electrolyzer"]
    storage = parameters.components["storage"]
    
    hourlymatching = (etc_length==1)
    println("Hourly matching: ", hourlymatching)
    maxcapgenerators = 100

    # Defining the variables
        # Design (these caps are currently arbitrary, should be added to json file)
    @variables(model, begin
        0 <= c_solar <= maxcapgenerators #Solar capacity
        0 <= c_wind_on <= maxcapgenerators #Onshore capacity
        0 <= c_wind_off <= maxcapgenerators #Offshore capacity
        0 <= c_electrolyzer <= maxcapelectrolyzer(parameters) #Electrolyzer capacity
        0 <= c_storage <= maxcapstorage(parameters) #Storage capacity
    end)

    # Optional constraints for fixed capacities
    # @constraints(model, begin
    #     c_solar == 37.4322
    #     c_wind_on == 15.0799
    #     c_wind_off == 0
    #     c_electrolyzer == 6.00945
    #     c_storage == 2918.96
    # end)

        # Operation
    @variables(model, begin
        0 <= p_DAM_sell[1:simulation_length] #Power sold to the DAM
        0 <= p_DAM_buy[1:simulation_length] # Power bought from the DAM
        0 <= p_electrolyzer[1:simulation_length] #Electrolyzer power
        os_electrolyzer[1:simulation_length], Bin #Electrolyzer on/standby status
        0 <= h_electrolyzer[1:simulation_length] #Electrolyzer hydrogen production
        0 <= h_storage_soc[1:simulation_length] #Electrolyzer hydrogen storage state of charge
        0 <= h_storage_in[1:simulation_length] #Electrolyzer hydrogen storage input
        0 <= h_storage_out[1:simulation_length] #Electrolyzer hydrogen storage output
        0 <= h_demand[1:simulation_length] #Hydrogen output
    end)

    if hourlymatching
        @constraint(
            model,
            [t in 1:simulation_length],
            p_DAM_buy[t] <= c_electrolyzer*electrolyzer.P_standby
        )
    # else
    #     @constraint(
    #         model,
    #         [t in 1:simulation_length],
    #         p_DAM_buy[t] <= 3*maxcapgenerators
    #     )
    end

    # #Putting bounds on buying and selling
    # @constraint(
    #     model,
    #     [t in 1:simulation_length],
    #     p_DAM_sell[t] <= 3*maxcapgenerators
    # )

        #Temporal correlation constraint
    @constraint(
        model, 
        [k in (1:etc_length:simulation_length)], 
        (sum((getnetpower(model, electrolyzer,t)) for t in (k:(k+etc_length-1))) <= 
        0)
    )
    
    # Defining the constraints
        # Number of upper bounds on variables that dramatically increase the model speed
            # Upper bound electrolyzer capacity, which should never exceed the max rated power of the RES
    @constraint(
        model, 
        (c_electrolyzer <= (c_wind_on+c_wind_off+c_solar))
    )

            # Upper bound electrolyzer hydrogen storage input, which should never exceed the max rated power of the electrolyzer
    @constraint(
        model,
        [t in 1:simulation_length],
        h_storage_in[t] <= hourly_target*5 #temporary fix, should be replaced by c_compressor or something of the sort
    )

    @constraint(
        model,
        [t in 1:simulation_length],
        h_storage_out[t] <= hourly_target*5
    )

        # Physical constraints
            # Energy balance
    @constraint(
        model, 
        [t in 1:simulation_length], 
        (p_electrolyzer[t] + p_DAM_sell[t] == 
        getenergyin(model, timeseries, t))
    )

        #  Upper bound electrolyzer power
    @constraint(
        model, 
        [t in 1:simulation_length], 
        (p_electrolyzer[t] <= 
        upperbound_electrolyzer(model, electrolyzer, t))
    ) 

        #  Lower bound electrolyzer power
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
    
        # Hydrogen storage initialization
    @constraint(
        model,
        (h_storage_soc[1] == storage.soc_init*c_storage + h_storage_in[1]*storage.efficiency - h_storage_out[1])
    )
    
        # Electrolyzer hydrogen storage state of charge
    @constraint(
        model,
        [t in 2:simulation_length],
        (h_storage_soc[t] == h_storage_soc[t-1] + h_storage_in[t]*storage.efficiency - h_storage_out[t])
    )

        # Hydrogen storage final state of charge equal to initial state of charge
    @constraint(
        model,
        h_storage_soc[1] == h_storage_soc[simulation_length]
    )
   
        # Hydrogen output
    @constraint(
        model,
        [t in 1:simulation_length],
        (h_demand[t] == h_electrolyzer[t] + h_storage_out[t] - h_storage_in[t])
    )

        # Hourly flexibility constraint
    @constraint(
        model,
        [t in 1:simulation_length],
        (1-(flexrange/(2-flexrange)))*hourly_target <=
        h_demand[t] <=
        (1+(flexrange/(2-flexrange)))*hourly_target 
    )

        # Hydrogen temporal correlation constraint
    @constraint(
        model, 
        [h in (1:htc_length:simulation_length)], 
        (sum(h_demand[t] for t in (h:(h+htc_length-1))) == 
        htc_length*hourly_target)
    )

    # Defining the objective function
    @objective(
        model, 
        Min, 
        (getinvestmentcosts(model, parameters) +
        sum(getmarketvalue(model, timeseries, lambda_TSO, t) for t in 1:simulation_length))
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

- `var_data::Dict{String,Any}`: Dictionary containing the results of the optimization model (dataframe and figures)
"""
function getresults(
    model::JuMP.Model, 
    parameters::Parameter_Data,
    verbose::Bool
)
    if termination_status(model) != MOI.OPTIMAL
        return NaN
    end

    var_data = getresultsdict(model)
    # Defining the most commonly used parameternames for readability

    capacitycomponentpairs = zip(
        (var_data[:c_solar], var_data[:c_wind_on], var_data[:c_wind_off], var_data[:c_electrolyzer], var_data[:c_storage]),
        (parameters.components[component] for component in ("solar", "wind_on", "wind_off", "electrolyzer", "storage"))
    )

    # Preparing the data for plotting
    p_solar = parameters.timeseries.solar*var_data[:c_solar]
    p_wind_on = parameters.timeseries.wind_on*var_data[:c_wind_on]
    p_wind_off = parameters.timeseries.wind_off*var_data[:c_wind_off]
    p_DAM = parameters.timeseries.price_DAM
    p_DAM_sell = var_data[:p_DAM_sell]

    # If the DAM price is not used, the DAM price is set to zero
    if haskey(var_data, :p_DAM_buy)
        p_DAM_buy = var_data[:p_DAM_buy]
    else
        p_DAM_buy = zeros(length(p_DAM))
    end

    power_data = DataFrame(
        :p_solar => p_solar,
        :p_wind_on => p_wind_on,
        :p_wind_off => p_wind_off,
        :p_electrolyzer => var_data[:p_electrolyzer],
        :p_DAM_buy => p_DAM_buy,
        :p_DAM_sell => p_DAM_sell
    )

    hydrogen_data = DataFrame(
        :h_electrolyzer => var_data[:h_electrolyzer],
        :h_storage_in => var_data[:h_storage_in],
        :h_storage_out => var_data[:h_storage_out],
        :h_demand => var_data[:h_demand],
        :h_storage_soc => var_data[:h_storage_soc]
    )

    # Extracting outcomes
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
        println(outcome_data)
    end

    return Dict(
        :power_data => power_data,
        :hydrogen_data => hydrogen_data,
        :outcome_data => outcome_data
    )
end

function savedata(
    data::Dict{Any,Any},
    parameters::Parameter_Data
)
    filename = "results_$(parameters.scenario.etc_length)_$(Parameters.scenario.htc_length)_$(parameters.scenario.flexrange).csv"
    # Saving the results to a csv file
    combined_data = hcat(data[:power_data], data[:hydrogen_data])
    CSV.write(filename, combined_data)
end

function addLCOH(
    results::DataFrame
)
    LCOH = []
    for (i, result) in enumerate(results[:, "Results"])
        if isa(result, Dict)
            append!(LCOH, result[:outcome_data][2, :cost]) 
        elseif isnan(result)
            append!(LCOH, NaN)
        end 
    end
    results[:, "LCOH"] = LCOH
    return results
end

function makeplot(
    flows::Array{Array{Float64,1},1},
    labels::Array{String,1};
    twindata::Union{Array{Float64,1}, Bool}=false,
    twinlabel::Union{String, Bool}=false,
    subsetsize::Union{Int64, Bool}=false
)
    # Removing the flows and labels that are all zeros
    zeros = findall(x->iszero(x), flows)
    deleteat!(flows, zeros)
    deleteat!(labels, zeros)

    # Taking flow subsets 
    if subsetsize != false
        println("Taking a subset of the flows")
        flows = map(x->x[1:subsetsize], flows)
        if isa(twindata, Array)
            twindata = twindata[1:subsetsize]
        end
    end

    fig = Figure()

    ax1 = Axis(fig[1,1])
    series!(ax1, flows, labels=labels)

    # Adding legend
    axislegend(ax1; position=:lt)

    # Adding a twin axis if needed
    if isa(twindata, Array)
        ax2 = Axis(fig[1,1], yaxisposition=:right)
        hidespines!(ax2)
        hidexdecorations!(ax2)
        series!(ax2, [twindata], labels=[twinlabel], solid_color=:black, linestyle=:dash)
        axislegend(ax2; position=:rt)
    end


    return fig, ax1
end

function powerfig(
    powerflow::Array{Array{Float64,1},1},
    powerlabels::Array{String,1};
    twindata::Union{Array{Float64,1}, Bool}=false,
    subsetsize::Union{Int64, Bool}=false
)
    fig, ax1 = makeplot(powerflow, powerlabels; twindata=twindata, subsetsize=subsetsize)
    ax1.title = "Power flows"
    ax1.ylabel ="Power [MW]"
    ax1.xlabel = "Time [h]"
    #ax2[:set_ylabel]("DAM price [€/MWh]")
    return fig
end

function hydrogenfig(
    hydrogenflow::Array{Array{Float64,1},1},
    hydrogenlabels::Array{String,1};
    twindata::Union{Array{Float64,1}, Bool}=false,
    subsetsize::Union{Int64, Bool}=false
)
    fig, ax1 = makeplot(hydrogenflow, hydrogenlabels; twindata=twindata, twinlabel="SOC", subsetsize=subsetsize)
    ax1.title = "hydrogen flows"
    ax1.ylabel ="hydrogen [kg]"
    ax1.xlabel = "Time [h]"
    #ax2[:set_ylabel]("DAM price [€/MWh]")
    return fig
end

function plotdata(
    data::Dict{Symbol,DataFrame}
)
    # set_theme!(;
    #     resolution=(800, 600),
    #     background_color="white",
    #     fontsize=8,
    #     Axis(;
    #         xgridstyle=:dash,
    #         ygridstyle=:dash,
    #         ),
    # )
    powerdata = data[:power_data]
    hydrogendata = data[:hydrogen_data]

    # Collecting all the power data
    powerlabels = ["Solar", "Wind (onshore)", "Wind (offshore)", "Buy", "Sell", "Electrolyzer"]
    powerflow = [powerdata[:, datapoint] for datapoint in [:p_solar, :p_wind_on, :p_wind_off, :p_DAM_buy, :p_DAM_sell, :p_electrolyzer]]

    # Collecting all the power data
    hydrogenlabels = ["Electrolyzer", "Storage in", "Storage out", "Demand"]
    hydrogenflow = [hydrogendata[:, datapoint] for datapoint in [:h_electrolyzer, :h_storage_in, :h_storage_out, :h_demand]]
    hsoc = hydrogendata[:, :h_storage_soc]

    # Plotting the power flows over the entire period
    println("Working on full plot")
    pfig = powerfig(powerflow, powerlabels)
    println("Working on subset plot")
    sub_pfig= powerfig(powerflow, powerlabels; subsetsize=(24*14))
    println("Working on hydrogen plot")
    hplot = hydrogenfig(hydrogenflow, hydrogenlabels)
    println("Working on hydrogen subplot")
    hsubplot = hydrogenfig(hydrogenflow, hydrogenlabels; twindata=hsoc, subsetsize=(24*14))
    
    return pfig, sub_pfig, hplot, hsubplot
end

function make3dplot(
    results::DataFrame
)
    X = convert(Array{Float32,1}, results[:, :etc_length])
    Y = convert(Array{Float32,1}, 100*results[:, :flexrange])
    Z = convert(Array{Float32,1}, results[:, :htc_length])

    C = [convert(Float32, x) for x in (results[:, :LCOH])]

    fig = Figure()
    ax1=Axis3(fig[1, 1], aspect=(1,1,1), azimuth = pi/4.8,elevation=pi/6)

    plt1 = scatter!(X, Y, Z, markersize=20, label="Scenario", color=C)
    #Legend(fig[1,2], ax1, valign=:top)
    Colorbar(fig[1,2], plt1, label="LCOH [€/kg]")
    
    ax1.xlabel="ETC length (h)"
    ax1.ylabel="Flexibility range (%)"
    ax1.zlabel="HTC length (h)"
    ax1.title="LCOH"
    fig
end

function makeheatmap(
    results::DataFrame;
    normalized::Bool=true
)
    X=:etc_length
    Y2=:flexrange
    Y1=:htc_length

    # Sorting results
    sort!(results, [X, Y1, Y2])

    timedict = Dict(
        1 => "Hourly",
        24 => "Daily",
        168 => "Weekly",
        336 => "Biweekly",
        672 => "Monthly",
        2016 => "Quarterly",
        4032 => "Half-yearly",
        8064 => "Yearly"
    )

    Xlabels = [timedict[label] for label in unique(results[:, X])]
    Y2labels = ["$(100*label)%" for label in unique(results[:, Y2])]
    Y1labels = [timedict[label] for label in unique(results[:, Y1])]
    Xscenarios = length(Xlabels)
    Yscenarios = length(Y1labels)*length(Y2labels)

    Xvals = repeat(range(1, Xscenarios), inner=convert(Int64,Yscenarios))
    Yvals = repeat(range(1, Yscenarios), convert(Int64,Xscenarios))

    Cvals = [convert(Float32, x) for x in (results[:, :LCOH])]
    if normalized
        Cvals = Cvals./mean(filter(!isnan, Cvals))
    end

    fig = Figure()
    ax1=Axis(
        fig[1, 1], 
        xlabel="ETC length", 
        ylabel="HTC length", 
        title="Average cost of hydrogen production for various scenarios without optimization [€/kg]",
        xticks=(range(1, Xscenarios),Xlabels),
        yticks=(
            range(
                start=length(Y2labels)/2+0.5, 
                stop=Yscenarios-length(Y2labels)/2+0.5, 
                length=length(Y1labels)
            ), 
            Y1labels
        ),
        xminorticks=range(start=0.5, stop=Xscenarios+0.5, length=Xscenarios+1),
        yminorticks=range(start=0.5, stop=Yscenarios+0.5, length=length(Y1labels)+1),
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminortickalign=0.75,
        yminortickalign=0.75,
        xminorticksize=60,
        yminorticksize=40,
        xticksize=0,
        yticksize=0
    )

    # Calculate colorbar range
    # Find maximum deviation from mean and round to nearest 0.1
    maxdeviation = ceil(maximum(abs.(Cvals.-1).*100))
    colorbarrange = (100-maxdeviation, 100+maxdeviation)./100
    colorbarticks = (100-maxdeviation:10:100+maxdeviation)./100
    colorbarticklabels = ["$(convert(Int64, round(100*tick)))%" for tick in colorbarticks]

    println(maxdeviation)
    println(colorbarrange)
    println([x for x in colorbarticks])
    println(colorbarticklabels)

    hm = heatmap!(Xvals, Yvals, Cvals, colormap=cgrad(:RdYlGn_8, rev=true), colorrange=colorbarrange, nan_color = :snow2)

    ax2=Axis(
        fig[1, 1],
        yaxisposition=:right,
        ylabel="Flexibility range (%)",
        yticks=(
            range(
                start=1/4,
                stop=10-1/4,
                length=Yscenarios
            ),
            repeat(Y2labels, length(Y1labels))
        ),
        ylabelsize=16,
        yticklabelsize=10,
        yminorticksvisible = true,
        yminortickalign=1,
        yminorticksize=5,
        yticksize=0,
        yticksmirrored=true,
        yminorgridvisible=true,
        yminorgridcolor=:black,
        yminorgridwidth=2,
    )
    hidespines!(ax2)
    hidexdecorations!(ax2)
    
    Colorbar(fig[:,2], hm, label="Deviations from mean [%]", ticks=(colorbarticks, colorbarticklabels))
fig
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
    optimizeptx(
        parameters::Parameter_Data, 
        verbose::Bool=true, 
        savepath::String="$(home_dir)/Results", 
        savefigs::Bool=true, 
        savecsv::Bool=true

This function caries out the optimization. This method works on the Parameter_Data struct.
    It takes a parameterfilename and solves the optimization model with the given parameters.

# Arguments

- `parameters::Parameter_Data`: Parameter data for the optimization model
- `verbose::Bool=true`: If true, prints the parameters to the console
- `savepath::String="$(home_dir)/Results"`: Path to save the results
- `savefigs::Bool=true`: If true, saves the figures to the savepath
- `savecsv::Bool=true`: If true, saves the csv files to the savepath

# Returns

- `results::Dict{String, Any}`: Dictionary containing the results of the optimization model

"""
function optimizeptx(
    parameters::Parameter_Data, 
    verbose::Bool=true,
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=true
)
    if verbose
        println("Parameters: ", parameters)
    end
    solvedmodel = solveoptimizationmodel(parameters, 120, 0.05)
    if termination_status(solvedmodel) == MOI.OPTIMAL
        println("Optimal solution found.")
        results = getresults(solvedmodel, parameters, verbose)
    else
        throw(ErrorException("Optimal solution not found."))
    end
    return results
end

"""
    optimizeptx(
        parameters::DataFrame, 
        verbose::Bool=true, 
        savepath::String="$(home_dir)/Results", 
        savefigs::Bool=true, 
        savecsv::Bool=true
    )

    Same as optimizeptx, but then the method for DataFrame as input. 
"""

function optimizeptx(
    parameterdf::DataFrame, 
    verbose::Bool=true,
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=true
)
    if verbose
        println("Parameters: ", parameterdf)
    end
    solutiondataframe = select(parameterdf, Not("Parameters"))
    solutiondataframe[:, "Results"] = [
        getresults(
            solveoptimizationmodel(
                parameters, 
                120, 
                0.05
            ), 
            parameters, 
            verbose
        ) for parameters in parameterdf[:, "Parameters"]
    ]
    return addLCOH(solutiondataframe)
end


function mainscript(
    parameterfilename::String,
    verbose::Bool=true,
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=true
)   
    parameters = fetchparameterdata(parameterfilename)

    results = optimizeptx(parameters, verbose, savepath)
    #plotdf(results)
    return results
end

