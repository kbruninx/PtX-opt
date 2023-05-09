
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
using StatsBase
using CairoMakie
using Tables
CairoMakie.activate!(type = "svg")

# Directory
const home_dir = @__DIR__


# Useful structures
struct Scenario_Data
    name::String
    # Financial parameters
    discount_rate::Float64
    lambda_TSO::Float64 # TSO penalty factor
    # Hydrogen target
    hourly_target::Float64 # kg/h of hydrogen production
    no_etc_periods::Int64 #number of hours per electricity temporal correlation period
    no_htc_periods::Int64 #number of hours per hyrdogen temporal correlation period
    flexrange::Float64 #flexibility range in % of the hourly target
    # Cap on combined generation capacity
    maxcapgeneration::Number
end

struct Component_Data
    name::String
    systemprice::Number    # €/unit
    fixopex::Number        # fixed operating expenses as a fraction of the system price
    lifetime::Int64         # lifetime in years
    maxcapacity::Number          # maximum capacity in MW
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

struct Compressor_Data
    # Storage data structure, containing additional storage data needed for the optimization model
    component_data::Component_Data
    # Additional storage parameters
    constant::Float64     # K constant for the compressor
end

#Forwards function calls to the component_data field for the Storage_Data structure
@forward((Compressor_Data, :component_data), Component_Data)

struct Timeseries_Data
    # Timeseries data structure, containing all the timeseries data needed for the optimization model
    solar::Array{Float64,1}
    wind_on::Array{Float64,1}
    wind_off::Array{Float64,1}
    price_DAM::Array{Float64,1}
    decision::Array{Float64,1}
    ordering::Array{Float64,2}
end

struct Parameter_Data
    # Parameter data structure, containing all the parameter data needed for the optimization model
    scenario::Scenario_Data
    components::Dict{String, Union{Component_Data, Electrolyzer_Data, Storage_Data, Compressor_Data}}
    timeseries::Timeseries_Data
end

#Printing methods
Base.show(io::IO, s::Scenario_Data) = print(io, "Scenario Data:
    Name = $(s.name),
    Discount rate = $(s.discount_rate),
    TSO tariff = $(s.lambda_TSO),
    Hourly production target = $(s.hourly_target),
    Electrical Temporal Correlation length = $(s.no_etc_periods),
    Hydrogen Temporal Correlation length = $(s.no_htc_periods),
    Flexibility range = $(s.flexrange),
    Max combined generation capacity = $(s.maxcapgeneration)"
)

Base.show(io::IO, c::Component_Data) = print(io, "Component Data:
    Name = $(c.name),
    System price = $(c.systemprice),
    Fixed OPEX = $(c.fixopex),
    Lifetime = $(c.lifetime),
    Max capacity = $(c.maxcapacity)"
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

Base.show(io::IO, s::Compressor_Data) = print(io, "Compressor Data:
    $(s.component_data),
    Constant = $(s.constant)"
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
function fetchparameterdata(
    parameterfile
    )
    d = JSON.parsefile(parameterfile)
    scenarios = combinescenarios(d["scenarios"])
    components = combinecomponents(d["sensitivity_analysis"], d["components"])
    timeseries = timeseriesfromdict(d["timeseriesfiles"])

    # If it is only a case study, return a Parameter_Data struct
    if scenarios isa Scenario_Data && components isa Dict
        return Parameter_Data(
                scenarios,
                components,
                timeseries
        )
    # If there are multiple scenarios, and no sensitivity analysis, return a DataFrame with the Parameter_Data structs
    elseif scenarios isa DataFrame && components isa Dict
        parameterdata = select(scenarios, Not("Scenario"))
        parameterdata[:, "Parameters"] = [Parameter_Data(
            scenario,
            components,
            timeseries
            ) for scenario in scenarios[:, "Scenario"]]
        return parameterdata
    # If there is a sensetivity analysis, and only one scenario, return a DataFrame with the Parameter_Data structs
    elseif scenarios isa Scenario_Data && components isa DataFrame
        parameterdata = select(components, Not("ComponentData"))
        parameterdata[:, "Parameters"] = [Parameter_Data(
            scenarios,
            case,
            timeseries
            ) for case in components[:, "ComponentData"]]
        return parameterdata
    else
        error("This type of data entry is currently not supported")
    end
end

function combinecomponents(
    sensitivity::Dict,
    components::Dict
)
    componentscenarios = DataFrame(
        "Component"=> String[], 
        "Variable" => String[], 
        "Adjustment" => Float64[],
        "ComponentData" => Dict{String, Any}[]
    )
    if sensitivity["variables"] isa Bool
        if sensitivity["variables"] == false
            return componentsfromdict(components)
        end
    end
    if sensitivity["variables"] isa Array
        for component in keys(components)
            for variable in sensitivity["variables"]
                for adjustment in sensitivity["adjustments"]
                    newcomponents = deepcopy(components)
                    newcomponents[component][variable] *= adjustment
                    push!(
                        componentscenarios,
                        (component, 
                            variable, 
                            adjustment, 
                            componentsfromdict(newcomponents)
                        ), 
                    )
                end
            end
        end
        return componentscenarios
    end
end

function permutate(
    dict::Dict{String, Any}
)
    return [Dict{String,Any}(zip(keys(dict), p)) for p in Iterators.product(values(dict)...)]
end

function permutatescenarios(
    dict::Dict{String, Any}
)
    #Checking for variables of interest (arrays)
    varsofinterest = [key for (key, value) in dict if value isa Array]
    #If there are none, return the scenario as a dict
    if varsofinterest == []
        return dict, varsofinterest
    end
    #Take out the name of scenario and permutate the rest
    name = pop!(dict, "name")
    permutations = permutate(dict)
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
    scenarios, varsofinterest = permutatescenarios(scenariodict)

    #If there are no variables of interest, return a Scenario_Data
    if varsofinterest == []
        return scenariofromdict(scenarios)
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
        dict["discount_rate"],
        dict["lambda_TSO"],
        dict["hourly_target"],
        dict["no_etc_periods"],
        dict["no_htc_periods"],
        dict["flexrange"],
        dict["maxcapgeneration"]
    )
end

function componentsfromdict(dict)
    out = Dict()
    for (k, v) in dict
        out[k] = Component_Data(
            k,
            v["systemprice"],
            v["fixopex"],
            v["lifetime"],
            v["maxcapacity"]
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
    out["compressor"] = Compressor_Data(
        out["compressor"],
        compressorconstant(dict["compressor"])/(1e3 * 3.6e3) # constant is kWs/kg -> MWh/kg
    )
    return out
end

# Timeseries manipulation
function timeseriesfromdict(
        timeseriesfiles::Dict{String, Any}
)
    timeseries_df = CSV.read(joinpath(home_dir, timeseriesfiles["profiles"]), DataFrame)
    decision = CSV.read(joinpath(home_dir, timeseriesfiles["decision"]), DataFrame)[:, "weights"]
    ordering = CSV.File(joinpath(home_dir, timeseriesfiles["ordering"])) |> Tables.matrix
    return Timeseries_Data(
        [timeseries_df[:, series] for series in ["Solar", "WindOnshore", "WindOffshore", "Price"]]...,
        decision,
        ordering
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
    parameters::Parameter_Data
)
    sum = 0
    c_compressor = model[:c_electrolyzer]*parameters.components["electrolyzer"].efficiency*parameters.components["compressor"].constant
    capacitycomponentpairs = zip(
        (model[char] for char in [:c_solar, :c_wind_on, :c_wind_off, :c_electrolyzer, :c_storage]),
        (parameters.components[component] for component in ("solar", "wind_on", "wind_off", "electrolyzer", "storage"))
    )
    for (capacity, component) in capacitycomponentpairs
        sum += getcomponentcosts(capacity, parameters.scenario, component)
        #sum += getcomponentcosts(c_compressor, parameters.scenario, parameters.components["compressor"])
    end
    return sum
end

function getcomponentcosts(
    capacity::Union{VariableRef, Float64}, 
    scenario::Scenario_Data,
    component::Union{Electrolyzer_Data, Storage_Data, Compressor_Data, Component_Data}
)
    crf = getcapitalrecoveryfactor(scenario.discount_rate, lifetime(component))
    return (
        systemprice(component)*capacity*
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
        max_storageoftotal*parameters.scenario.hourly_target*size(parameters.timeseries.ordering, 1)
    )
end

function compressorconstant(
    cd::Dict #Compressor dictionary
)
    R = 8.31446261815324 #J/(mol K)
    gamma = 1.4
    K = (
        ((R*cd["inlet_temperature"])/(2*(gamma-1)*cd["efficiency"]))*
        ((cd["outlet_pressure"]/cd["inlet_pressure"])^((gamma-1)/gamma) - 1)
    )
    return K
end

# Transforms an array of representative values to an array of hourly values over the year based on an ordering matrix
function rep_to_ext(
    repvals,
    ordering::Array{Float64, 2};
    daily::Bool = false
    
)
    #When the resolution of the sums is daily, we can speed up the calculation by first summing the representative values per day
    if daily == true
        reparray = [sum(repvals[i:i+23]) for i in 1:24:length(repvals)]
    #Otherwise we are going to reshape the representative values to a matrix of 24 columns and as many rows as there are representative days
    else 
        reparray = reshape(repvals, 24, length(repvals)÷24)'
    end
    #We then perform a matrix multiplication which will yield a matrix with 24 columns and as many rows as there are total days, which we flatten to a vector
    return vec((ordering*reparray)')
end

# Generates a vector of vectors each representing a period out of a number of periods over a number of total days
# In the case of larger than daily periodes, we have to ensure the periods are sliced per day rather than per hour
function hourlyperiodsarray(
    periods::Int64;
    daysintotal::Int64=365
)
    # Checks if we're slicing days or hours
    if periods >= daysintotal
        # Checks if we can slice days evenly
        if (daysintotal*24)%periods != 0
            throw(ArgumentError("Can not slice days unevenly for $periods periods: $((daysintotal*24)%periods)"))
        end
        # Finds the indices of the start of each of the periods for an hourly resolution
        hourlyidx = [round(Int64, i) for i in range(0,daysintotal*24,periods+1)]
        # Fills out the idx array, mapping each start index to an array of hours for that entire period
        hourlyarray = [collect(hourlyidx[i-1]+1:hourlyidx[i]) for i in 2:length(hourlyidx)]
        return hourlyarray
    else
        # Finds the indices of the start of each of the periods for a daily resolution
        dailyidx = [round(Int64, i) for i in range(0,daysintotal,periods+1)]
        # Fills out the idx array, mapping each start index to an array of days for that entire period
        dailyarray = [collect(dailyidx[i-1]:dailyidx[i]-1) for i in 2:length(dailyidx)]
        # Maps each day to an array of hours for that day
        return [vcat([collect(24*day+1:24*day+24) for day in days]...) for days in dailyarray]
    end
end

"""
    solveoptimizationmodel(
        parameters::Parameter_Data,
        time_limit::Int64,
        mip_gap::Float64
    )

This function takes parameter data and returns a solved model with the results for the given parameters.

# Arguments

- `parameters::Parameter_Data`: Structure containing all the parameters for the optimization model
- `time_limit::Int64`: Time limit in seconds for the optimization model
- `mip_gap::Float64`: MIP gap for the optimization model
    
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
    scenario = parameters.scenario # Scenario data object
    hourly_target = scenario.hourly_target # Hourly target of hydrogen production in kg
    no_etc_periods = scenario.no_etc_periods # Number of periods for the electrical temporal correlation
    no_htc_periods = scenario.no_htc_periods # Number of periods for the hydrogen temporal correlation
    flexrange = scenario.flexrange # Flexibility as a percentage deviation from the maximum downstream hydrogen inlet
    lambda_TSO = scenario.lambda_TSO # Tariff for energy purchased from the grid

    components = parameters.components
    electrolyzer = components["electrolyzer"]
    storage = components["storage"]
    compressor = components["compressor"]
    
    timeseries = parameters.timeseries # Timeseries data object
    decision_array = repeat(timeseries.decision, inner=24) # Array of weights for the representative days
    representative_timesteps = size(timeseries.ordering, 2)*24 # Number of timesteps in the representative days
    total_days = size(timeseries.ordering, 1) # Number of days in the total timeseries
    total_timesteps = total_days*24  # Number of timesteps in the total timeseries

    println("Putting together the arrays with the periods: $no_htc_periods and $total_days")
    htc_periods = hourlyperiodsarray(no_htc_periods, daysintotal=total_days) # Array with the periods for the hydrogen temporal correlation
    etc_periods = hourlyperiodsarray(no_etc_periods, daysintotal=total_days) # Array with the periods for the electrical temporal correlation

    hourlymatching = (no_etc_periods==8760) # Boolean for hourly matching, speed ups are implemented if true
    println("Hourly matching: ", hourlymatching)

    # Defining the variables
        # Design (these caps are currently arbitrary, should be added to json file)
    @variables(model, begin
        0 <= c_solar #Solar capacity
        0 <= c_wind_on #Onshore capacity
        0 <= c_wind_off #Offshore capacity
        0 <= c_electrolyzer #<= maxcapelectrolyzer(parameters) #Electrolyzer capacity
        0 <= c_compressor #Compressor capacity
        0 <= c_storage #<= maxcapstorage(parameters) #Storage capacity
    end)

    # Limiting capacities to a combined maximum
    # First the combined generation capacity
    if !isa(scenario.maxcapgeneration, Bool)
        @constraint(model, c_solar + c_wind_on + c_wind_off <= scenario.maxcapgeneration)
    end

    # Then each individual generation capacity, if given
    for (string, symbol) in zip(("solar", "wind_on", "wind_off"), (:c_solar, :c_wind_on, :c_wind_off))
        if !isa(components[string].maxcapacity, Bool)
            @constraint(model, model[symbol] <= components[string].maxcapacity)
        end
    end

    # Optional constraints for fixed capacities
    # @constraints(model, begin
    #     c_solar == 10
    #     c_wind_on == 0
    #     c_wind_off == 30
    #     c_electrolyzer == 7
    # end)

        # Operation
    @variables(model, begin
        0 <= p_DAM_sell[1:representative_timesteps] #Power sold to the DAM in MW
        0 <= p_DAM_buy[1:representative_timesteps] # Power bought from the DAM in MW
        0 <= p_electrolyzer[1:representative_timesteps] #Electrolyzer power in MW
        0 <= p_compressor[1:representative_timesteps] #Compressor power in MW
        os_electrolyzer[1:representative_timesteps], Bin #Electrolyzer on/standby status
        0 <= h_electrolyzer[1:representative_timesteps] #Electrolyzer hydrogen production in kg
        0 <= h_storage_soc[1:representative_timesteps] #Electrolyzer hydrogen storage state of charge in kg
        0 <= h_storage_in[1:representative_timesteps] #Electrolyzer hydrogen storage input in kg
        0 <= h_storage_out[1:representative_timesteps] #Electrolyzer hydrogen storage output in kg
        0 <= h_demand_dir[1:representative_timesteps] #Hydrogen output directly in kg
        0 <= h_demand[1:representative_timesteps] #Hydrogen output in kg
    end)

    # Defining the constraints

    # If hourly matching, we limit the potential power bought from the DAM s.t. we can't use it for electrolysis
    if hourlymatching
        @constraint(
            model,
            [t in 1:representative_timesteps],
            p_DAM_buy[t] <= c_electrolyzer*electrolyzer.P_standby+p_compressor[t]
        )
    end

    # Defining the net consumption as the power bought, minus all the power sold and used for non-electrolyis purposes
    netconsumption_rep = [
        p_DAM_buy[t] - (
            p_DAM_sell[t] +
            os_electrolyzer[t]*p_electrolyzer[t]*electrolyzer.efficiency*compressor.constant +
            (1 - os_electrolyzer[t]) * (electrolyzer.P_standby * c_electrolyzer))
        for t in 1:representative_timesteps
    ]

    # Extending the net consumption to the total timeseries
    netconsumption_ext = rep_to_ext(netconsumption_rep, timeseries.ordering)

        #Temporal correlation constraint
    @constraint(
        model, 
        [period in etc_periods], 
        (sum(
             netconsumption_ext[t]
        for t in period) <= 0)
    )
    
        # Number of upper bounds on variables that increase the model speed
            # Upper bound electrolyzer capacity, which should never exceed the max rated power of the RES
    @constraint(
        model, 
        (c_electrolyzer <= (c_wind_on+c_wind_off+c_solar))
    )

        # Physical constraints
            # Energy balance
    @constraint(
        model, 
        [t in 1:representative_timesteps], 
        (p_electrolyzer[t] + p_compressor[t] + p_DAM_sell[t] == (
        c_solar*timeseries.solar[t] +
        c_wind_on*timeseries.wind_on[t] +
        c_wind_off*timeseries.wind_off[t] + 
        p_DAM_buy[t]))
    )

        #  Upper bound electrolyzer power
    @constraint(
        model, 
        [t in 1:representative_timesteps], 
        (p_electrolyzer[t] <= (
        c_electrolyzer*os_electrolyzer[t] +
        electrolyzer.P_standby*c_electrolyzer*(1-os_electrolyzer[t])))
    )

        #  Lower bound electrolyzer power
    @constraint(
        model, 
        [t in 1:representative_timesteps], 
        (p_electrolyzer[t] >= (
        c_electrolyzer*electrolyzer.P_min*os_electrolyzer[t] +
        electrolyzer.P_standby*c_electrolyzer*(1-os_electrolyzer[t])))
    )

        # Electrolyzer hydrogen production
    @constraint(
        model,
        [t in 1:representative_timesteps],
        (h_electrolyzer[t] == (
            electrolyzer.efficiency * 
            p_electrolyzer[t] * 
            os_electrolyzer[t]))
    )

    h_storage_soc_ext = rep_to_ext([h_storage_soc[t] for t in 1:representative_timesteps], timeseries.ordering)
    h_storage_in_ext = rep_to_ext([h_storage_in[t] for t in 1:representative_timesteps], timeseries.ordering)
    h_storage_out_ext = rep_to_ext([h_storage_out[t] for t in 1:representative_timesteps], timeseries.ordering)
    h_demand_ext = rep_to_ext([h_demand[t] for t in 1:representative_timesteps], timeseries.ordering)

        # Hydrogen storage capacity
    @constraint(
        model,
        [t in 1:representative_timesteps],
        (h_storage_soc[t] <= c_storage)
    )
    
        # Hydrogen storage initialization
    @constraint(
        model,
        (h_storage_soc_ext[1] == storage.soc_init*c_storage + h_storage_in_ext[1] - h_storage_out_ext[1])
    )
    
        # Electrolyzer hydrogen storage state of charge
    @constraint(
        model,
        [t in 2:total_timesteps],
        (h_storage_soc_ext[t] == h_storage_soc_ext[t-1] + h_storage_in_ext[t] - h_storage_out_ext[t])
    )

        # Hydrogen storage final state of charge equal to initial state of charge
    @constraint(
        model,
        h_storage_soc_ext[1] == h_storage_soc_ext[total_timesteps]
    )

    # Compressor power consumption
    @constraint(
        model,
        [t in 1:representative_timesteps],
        p_compressor[t] == compressor.constant * h_storage_in[t]
    )

    @constraint(
        model,
        [t in 1:representative_timesteps],
        p_compressor[t] <= c_compressor
    )
   
    # Electrolyzer output
    @constraint(
        model,
        [t in 1:representative_timesteps],
        (h_electrolyzer[t] == h_demand_dir[t] + h_storage_in[t])
    )

    # System output
    @constraint(
        model,
        [t in 1:representative_timesteps],
        (h_demand[t] == h_demand_dir[t] + h_storage_out[t])
    )

    #    Hourly flexibility constraint
    @constraint(
        model,
        [t in 1:representative_timesteps],
        (1-(flexrange/(2-flexrange)))*hourly_target <=
        h_demand[t] <=
        (1+(flexrange/(2-flexrange)))*hourly_target 
    )

    # Hydrogen temporal correlation constraint
    @constraint(
        model, 
        [period in htc_periods],
        (sum(h_demand_ext[t] for t in period) == 
        length(period)*hourly_target)
    )

    # Defining the objective function
    @objective(
        model, 
        Min, 
        (
            getinvestmentcosts(model, parameters) +
            sum(
                decision_array[t] * (
                p_DAM_buy[t] * (timeseries.price_DAM[t] + lambda_TSO) - 
                p_DAM_sell[t] * timeseries.price_DAM[t])
                for t in 1:representative_timesteps)
        )
    )

    # Solving the model
    optimize!(model)

    return model
end

function unconstrainedmodel(
    parameters::Parameter_Data, 
    time_limit::Int64
)
    # Defining the model
    model = Model(Gurobi.Optimizer)

    # Setting the time limit and the mip gap
    set_time_limit_sec(model, time_limit)
    # set_optimizer_attribute(model, "DualReductions", 0)

    # Defining some variables for visibility
    scenario = parameters.scenario
    components = parameters.components
    timeseries = parameters.timeseries

    decision_array = repeat(timeseries.decision, inner=24)
    representative_timesteps = size(timeseries.ordering, 2)*24

    @variables(model, begin
        0 <= c_solar #Capacity of solar power plant
        0 <= c_wind_on #Capacity of wind power plant
        0 <= c_wind_off #Capacity of wind power plant
    end)

    # Limitting max capacities
    if !isa(scenario.maxcapgeneration, Bool)
        @constraint(model, c_solar + c_wind_on + c_wind_off == scenario.maxcapgeneration)
    end

    # Limiting individual capacities
    for (string, symbol) in zip(("solar", "wind_on", "wind_off"), (:c_solar, :c_wind_on, :c_wind_off))
        if !isa(components[string].maxcapacity, Bool)
            @constraint(model, model[symbol] <= components[string].maxcapacity)
        end
    end
    
    capacitycomponentpairs = zip(
        (c_solar, c_wind_on, c_wind_off),
        (components["solar"], components["wind_on"], components["wind_off"])
    )

    @objective(
        model,
        Min,
        sum(getcomponentcosts(capacity, scenario, component) for (capacity, component) in capacitycomponentpairs) - 
        sum(decision_array[t] * timeseries.price_DAM[t] * energygenerated(model, timeseries, t) for t in 1:representative_timesteps)
    )

    optimize!(model)

    return model

end

function getresultsdict(m::JuMP.Model)
    return Dict(
        k => value.(v) for 
        (k, v) in object_dictionary(m)
    )
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
    parameters::Parameter_Data;
    verbose::Bool,
    opportunity_cost::Union{Float64, Int} = 0
)
    if termination_status(model) == MOI.INFEASIBLE
        return "Infeasible"
    elseif termination_status(model) == MOI.TIME_LIMIT
        return "Time limit reached"
    elseif termination_status(model) != MOI.OPTIMAL
        return "No optimal solution found for unknown reason"
    end

    # Defining the most commonly used parameternames for readability
    var_data = getresultsdict(model)
    total_timesteps = size(parameters.timeseries.ordering, 1)*24
    decision_array = repeat(parameters.timeseries.decision, inner=24)
    ordering = parameters.timeseries.ordering

    capacitycomponentpairs = zip(
        (var_data[char] for char in (:c_solar, :c_wind_on, :c_wind_off, :c_electrolyzer, :c_storage, :c_compressor)),
        (parameters.components[component] for component in ("solar", "wind_on", "wind_off", "electrolyzer", "storage", "compressor"))
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

    total_cost = objective_value(model)-opportunity_cost
    total_profit_DAM = sum(decision_array.*(p_DAM_sell.*p_DAM.-(p_DAM_buy.*(p_DAM.+parameters.scenario.lambda_TSO))))
    capital_cost = objective_value(model)+total_profit_DAM

    # Extracting outcomes
    outcome_data = Dict{String, Union{Float64, Dict{String, Float64}}}(
        "Total costs" => total_cost,
        "Unadjusted costs" => objective_value(model),
        "Average cost of Hydrogen" => (ototal_cost)/(parameters.scenario.hourly_target*total_timesteps),
        "Unadjusted average cost of Hydrogen" => objective_value(model)/(parameters.scenario.hourly_target*total_timesteps),
        "Electrolyzer capacity factor" => mean((var_data[:p_electrolyzer]./var_data[:c_electrolyzer]), weights(decision_array)),
        "Total profit power market" => total_profit_DAM 
    )

    for (capacity, component) in capacitycomponentpairs
        componentcost = getcomponentcosts(capacity, parameters.scenario, component)
        outcome_data["$(name(component))"] = Dict(
            "capacity"=> capacity, 
            "cost" => componentcost, 
            "percentage" => componentcost/capital_cost
        )
    end


    # Printing the results if verbose
    if verbose
        print(json(outcome_data,4))
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
    filename = "results_$(parameters.scenario.no_etc_periods)_$(Parameters.scenario.no_htc_periods)_$(parameters.scenario.flexrange).csv"
    # Saving the results to a csv file
    combined_data = hcat(data[:power_data], data[:hydrogen_data])
    CSV.write(filename, combined_data)
end

function addoutcometodf(
    results::DataFrame,
    outcomes::Vector{Union{String}}
)
    for outcome in outcomes
        println("Adding outcome: $(outcome)")
        data=[]
        for result in results[:, "Results"]
            if result isa Dict{Symbol, Any}
                push!(data, result[:outcome_data][outcome])
            elseif  !(result isa Number) || isnan(result) 
                push!(data, NaN)
            end
        end
        dictdata = filter(t -> (typeof(t)==Dict{String, Float64}), data)
        if !isempty(dictdata)
            for key in keys(dictdata[1])
                insertcols!(results, "$(outcome)_$(key)" => [if typeof(datapoint)==Dict{String, Float64} datapoint[key] else datapoint end for datapoint in data])
            end
        else
            insertcols!(results, "$(outcome)" => data)
        end
    end
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
        start = 1#round(Int, length(flows[1])/2)
        flows = map(x->x[start:start+subsetsize], flows)
        if isa(twindata, Array)
            twindata = twindata[1:subsetsize]
        end
    end

    fig = Figure()

    ax1 = Axis(fig[1,1])
    series!(ax1, flows, labels=labels)

    # Adding legend
    axislegend(ax1; position=:lt)

    ax2 = Axis(fig[1,1], yaxisposition=:right)
    hidespines!(ax2)
    hidexdecorations!(ax2)

    # Adding a twin axis if needed
    if isa(twindata, Array)
        series!(ax2, [twindata], labels=[twinlabel], solid_color=:black, linestyle=:dash)
        axislegend(ax2; position=:rt)
        return fig, ax1, ax2
    end


    return fig, ax1, ax2
end

function powerfig(
    powerflow::Array{Array{Float64,1},1},
    powerlabels::Array{String,1};
    twindata::Union{Array{Float64,1}, Bool}=false,
    subsetsize::Union{Int64, Bool}=false
)
    fig, ax1, ax2 = makeplot(powerflow, powerlabels; twindata=twindata, twinlabel="DAM Price", subsetsize=subsetsize)
    ax1.title = "Power flows"
    ax1.ylabel ="Power [MW]"
    ax1.xlabel = "Time [h]"
    ax2.ylabel = "DAM price [€/MWh]"
    return fig
end

function hydrogenfig(
    hydrogenflow::Array{Array{Float64,1},1},
    hydrogenlabels::Array{String,1};
    twindata::Union{Array{Float64,1}, Bool}=false,
    subsetsize::Union{Int64, Bool}=false
)
    fig, ax1, ax2 = makeplot(hydrogenflow, hydrogenlabels; twindata=twindata, twinlabel="SOC", subsetsize=subsetsize)
    ax1.title = "Hydrogen flows"
    ax1.ylabel ="Hydrogen massflow [kg]"
    ax1.xlabel = "Time [h]"
    ax2.ylabel = "Hydrogen stored [kg]"
    return fig
end

function plotdata(
    data::Dict{Symbol,Any};
    parameters=false
)
    powerdata = data[:power_data]
    hydrogendata = data[:hydrogen_data]

    # Collecting all the power data
    powerlabels = ["Solar", "Wind (onshore)", "Wind (offshore)", "Buy", "Sell", "Electrolyzer"]
    powerflow = [powerdata[:, datapoint] for datapoint in [:p_solar, :p_wind_on, :p_wind_off, :p_DAM_buy, :p_DAM_sell, :p_electrolyzer]]
    if !(parameters isa Bool)
        damprices = powerdata[:, :p_DAM_sell]
    else 
        damprices = zeros(length(powerdata[:, :p_solar]))
    end

    # Collecting all the hydrogen data
    hydrogenlabels = ["Electrolyzer", "Storage in", "Storage out", "Demand"]
    hydrogenflow = [hydrogendata[:, datapoint] for datapoint in [:h_electrolyzer, :h_storage_in, :h_storage_out, :h_demand]]
    hsoc = hydrogendata[:, :h_storage_soc]

    # Plotting the flows
    println("Working on full plot")
    pfig = powerfig(powerflow, powerlabels)
    println("Working on subset plot")
    sub_pfig= powerfig(powerflow, powerlabels; subsetsize=(24*14))
    println("Working on hydrogen plot")
    hplot = hydrogenfig(hydrogenflow, hydrogenlabels; twindata=hsoc)
    println("Working on hydrogen subplot")
    hsubplot = hydrogenfig(hydrogenflow, hydrogenlabels; twindata=hsoc, subsetsize=(24*14))
    
    return pfig, sub_pfig, hplot, hsubplot
end

function make3dplot(
    results::DataFrame
)
    X = convert(Array{Float32,1}, results[:, :no_etc_periods])
    Y = convert(Array{Float32,1}, 100*results[:, :flexrange])
    Z = convert(Array{Float32,1}, results[:, :no_htc_periods])

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

#Unfinished function
function plotsensitivity(
    results::DataFrame;
    normalized::Bool=true,
    LCXY::Array{Char, 1}=["Variable", "Component", "Adjustment", "LCOH"]
)
    L = unique(results[LCXY[1]])
    C = []
    C = [[component[LCXY[4]]]]
    X=sort(unique(results[X]))

end

"""
    plotheatmap(
        results::DataFrame;
        C::Union{Symbol,String}=:LCOH,
        Clabel::String="LCOH [€/kg]",
        normalized::Bool=false,
        percentage::Bool=false
    )

Plots a heatmap of the results of a scenario analysis. 
The heatmap is plotted using the `C` column of the `results` DataFrame. 
The `Clabel` is the label of the colorbar. 
The `normalized` argument determines whether the heatmap is normalized. 
The `percentage` argument determines whether the heatmap is plotted in percentages.
"""
function plotheatmap(
    results::DataFrame;
    C::Union{Symbol,String}=:LCOH,
    Clabel::String="LCOH [€/kg]",
    normalized::Bool=false,
    percentage::Bool=false
)
    X=:no_etc_periods
    Y2=:flexrange
    Y1=:no_htc_periods

    # Sorting results
    sort!(results, [X, Y1, Y2])

    # Dictionary to map numbers to labels
    timedict = Dict(
        1 => "Yearly",
        4 => "Quarterly",
        12 => "Monthly",
        26 => "Biweekly",
        52 => "Weekly",
        365 => "Daily",
        8760 => "Hourly"
    )

    # Creating the labels
    Xlabels = [timedict[label] for label in unique(results[:, X])]
    Y2labels = ["$(100*label)%" for label in unique(results[:, Y2])]
    Y1labels = [timedict[label] for label in unique(results[:, Y1])]
    Xscenarios = length(Xlabels)
    Yscenarios = length(Y1labels)*length(Y2labels)

    # Creating the value arrays
    Xvals = repeat(range(1, Xscenarios), inner=convert(Int64,Yscenarios))
    Yvals = repeat(range(1, Yscenarios), convert(Int64,Xscenarios))
    Cvals = [convert(Float32, x) for x in (results[:, C])]

    # Normalizing the values if required
    if percentage
        Cvals = [Cvals[i]/results[:total_cost]]
    end
    if normalized
        Cvals = Cvals./mean(filter(!isnan, Cvals))
        Clabel = "Deviation from the mean [%]"
    end

    fig = Figure()
    ax1=Axis(
        fig[1, 1], 
        xlabel="ETC length", 
        ylabel="HTC length", 
        title="Average cost of hydrogen production for various scenarios with optimization [€/kg]",
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

    hm = heatmap!(Xvals, Yvals, Cvals, colormap=cgrad(:RdYlGn_8, rev=true), nan_color = :snow2)

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
    
    Colorbar(fig[:,2], hm, label=Clabel)
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
    parameters::Parameter_Data;
    verbose::Bool=true,
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=true,
    opportunity::Bool=true
)
    if verbose
        println("Parameters: ", parameters)
    end
    solvedmodel = solveoptimizationmodel(parameters, 600, 0.05)
    println("Opportunity is $(opportunity)")
    if opportunity
        solvedopportunity = unconstrainedmodel(parameters, 60)
        opportunity_cost = objective_value(solvedopportunity)
        println("Opportunity cost: ", opportunity_cost)
    else
        opportunity_cost = 0.0
    end
    if termination_status(solvedmodel) == MOI.OPTIMAL
        println("Optimal solution found.")
        results = getresults(
            solvedmodel, 
            parameters, 
            verbose=verbose, 
            opportunity_cost=opportunity_cost)
    else
        throw(ErrorException("Optimal solution not found. \n Termination status: $(termination_status(solvedmodel))"))
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
    parameterdf::DataFrame;
    verbose::Bool=true,
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=true,
    opportunity::Bool=true
)
    if verbose
        println("Parameters: ", parameterdf)
    end
    solutiondataframe = select(parameterdf, Not("Parameters"))
    results = []
    for parameters in parameterdf[:, "Parameters"]
        solvedmodel = solveoptimizationmodel(parameters, 600, 0.05)
        if opportunity
            opportunity_cost = objective_value(unconstrainedmodel(parameters, 60))
        else
            opportunity_cost = 0.0
        end
        push!(
            results,
            getresults(
                solvedmodel, 
                parameters, 
                verbose=verbose,
                opportunity_cost=opportunity_cost)
        )
    end
    solutiondataframe[:, "Results"] = results
    return solutiondataframe
end


function mainscript(
    parameterfilename::String,
    verbose::Bool=true;
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=false,
    opportunity::Bool=true
)   
    parameters = fetchparameterdata(parameterfilename)
    results = optimizeptx(parameters, verbose=verbose, savepath=savepath, savecsv=savecsv, opportunity=opportunity)
    return results
end

