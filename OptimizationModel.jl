
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
using JuMP, Gurobi, Clp # Optimization packages
using DataFrames, LinearAlgebra, Tables # Data packages
using Dates, ReusePatterns, Setfield # Utility packages
using Random, Distributions, StatsBase, Statistics # Statistics packages
using FileIO, JSON, CSV, JLD2 # File handling packages
using CairoMakie # Plotting packages
CairoMakie.activate!(type = "png")

# Directory
const home_dir = @__DIR__

year = 2020

# Useful structures
struct Scenario_Data
    name::String
    # Financial parameters
    discount_rate::Float64
    lambda_TSO::Float64 # TSO penalty factor
    # Hydrogen target
    hourly_target::Float64 # kg/h of hydrogen production
    no_etc_periods::Union{Int64, Nothing} #number of hours per electricity temporal correlation period
    no_htc_periods::Int64 #number of hours per hyrdogen temporal correlation period
    MPLDSP::Float64 #Minimum partial load of the DSP
    CFDSP::Float64 #Capacity factor of the DSP
    EI_limit::Union{Float64, Nothing} #Upper bound of the emission intensity in gCO2/kgH2
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
    solar::Vector{Float64}
    wind_on::Vector{Float64}
    wind_off::Vector{Float64}
    price_DAM::Vector{Float64}
    EI::Vector{Float64}
    representativedaysarray::Vector{Int64}
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
    Minimum Partial Load DSP = $(s.MPLDSP),
    Minimum Capacity factor DSP = $(s.CFDSP),
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
    Solar capacity factor mean = $(mean(t.solar, weights(repeat(decisionarray(t.ordering), inner=24))))
    Onshore wind capacity factor mean =  $(mean(t.wind_on, weights(repeat(decisionarray(t.ordering), inner=24)))),
    Offshore wind capacity factor mean =  $(mean(t.wind_off, weights(repeat(decisionarray(t.ordering), inner=24)))),
    Price DAM mean= $(mean(t.price_DAM, weights(repeat(decisionarray(t.ordering), inner=24)))),
    Representative days = $(size(t.ordering,1)) == $(length(decisionarray(t.ordering))) == $(round(Int64, length(t.price_DAM)/24)),
    Total days = $(size(t.ordering,2))"
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
        parameterdata = select(scenarios, Not("scenario"))
        parameterdata[:, "parameters"] = [Parameter_Data(
            scenario,
            components,
            timeseries
            ) for scenario in scenarios[:, "scenario"]]
        return parameterdata
    # If there is a sensetivity analysis, and only one scenario, return a DataFrame with the Parameter_Data structs
    elseif scenarios isa Scenario_Data && components isa DataFrame
        parameterdata = select(components, Not("component_data"))
        parameterdata[:, "parameters"] = [Parameter_Data(
            scenarios,
            case,
            timeseries
            ) for case in components[:, "component_data"]]
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
        "component"=> String[], 
        "variable" => String[], 
        "adjustment" => Float64[],
        "component_data" => Dict{String, Any}[]
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
    # 
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

    notvars = [(key, value) for (key,value) in dict if isnothing(value)]
    for varvalue in notvars
        delete!(dict, varvalue[1])
    end
    permutations = permutate(dict)
    #Put the name back in
    for p in permutations
        p["name"] = name
        for varvalue in notvars
            p[varvalue[1]] = varvalue[2]
        end
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
    df = DataFrame([columnname => [] for columnname in [varsofinterest..., "scenario"]])
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
        dict["MPLDSP"],
        dict["CFDSP"],
        dict["EI_limit"],
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
    ordering = transpose(CSV.File(joinpath(home_dir, timeseriesfiles["ordering"])) |> Tables.matrix)
    return Timeseries_Data(
        [timeseries_df[:, series] for series in ["Solar", "WindOnshore", "WindOffshore", "Price", "EI"]]...,
        unique(timeseries_df[:, "period"]),
        ordering
        )
end

function datetimefromrepdays(
    repdays::Vector{Int64}
)
    return [(DateTime(year, 1, 1, hour) +Dates.Day(day-1)) for day in repdays for hour in 0:23]
end

function datetimewholeyear(
    year::Int64=2020,
    numberofdays::Int64=365
)
    return [(DateTime(year, 1, 1, hour) +Dates.Day(day-1)) for day in 1:numberofdays for hour in 0:23]
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
    parameters::Parameter_Data;
    except = []
)
    sum = 0
    capacities = [:c_solar, :c_wind_on, :c_wind_off, :c_electrolyzer, :c_storage, :c_compressor]
    componentstr = ["solar", "wind_on", "wind_off", "electrolyzer", "storage", "compressor"]

    #Getting rid of exceptions
    if !isempty(except)
        exidx = findall(x -> x in except, componentstr)
        deleteat!(capacities, exidx)
        deleteat!(componentstr, exidx)
    end
    
    capacitycomponentpairs = zip(
        (model[char] for char in capacities),
        (parameters.components[component] for component in componentstr)
    )
    for (capacity, component) in capacitycomponentpairs
        sum += getcomponentcosts(capacity, parameters.scenario, component)
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

function ext_to_rep(
    extvals,
    ordering;

)
    #We first reshape the extended values to a matrix of 24 columns and as many rows as there are total days
    extarray = reshape(extvals, 24, length(extvals)÷24)
    #We then perform a matrix division which will yield a matrix with 24 columns and as many rows as there are representative days, 
    reparray = extarray/ordering
    # In order to flatten to a vector, we need to transpose the matrix first, because Julia flattens columnwise
    out=vec(reparray)    
    display(out)
    return out
end

# Transforms an array of representative values to an array of hourly values over the year based on an ordering matrix
function rep_to_ext(
    repvals,
    ordering::AbstractMatrix;
    daily::Bool = false
)
    #When the resolution of the sums is daily, we can speed up the calculation by first summing the representative values per day
    if daily == true
        reparray = permutedims([sum(repvals[i:i+23]) for i in 1:24:length(repvals)])
    #Otherwise we are going to reshape the representative values to a matrix of 24 columns and as many rows as there are representative days
    else 
        reparray = reshape(repvals, 24, length(repvals)÷24)
    end
    #We then perform a matrix multiplication which will yield a matrix with 24 columns and as many rows as there are total days, 
    # In order to flatten to a vector, we need to transpose the matrix first, because Julia flattens columnwise
    return vec(reparray*ordering)
end

# Generates a vector of vectors each representing a period out of a number of periods over a number of total days
# In the case of larger than daily periodes, we have to ensure the periods are sliced per day rather than per hour
function periodsarray(
    periods::Int64;
    daysintotal::Int64=365,
    daily::Bool=false
)
    # Checks if we're slicing days or hours
    if periods > daysintotal
        # Checks if we can slice days evenly
        if (daysintotal*24)%periods != 0
            throw(ArgumentError("Can not slice $daysintotal days unevenly for $periods periods: $((daysintotal*24)%periods)"))
        end
        # Finds the indices of the start of each of the periods for an hourly resolution
        hourlyidx = [round(Int64, i) for i in range(0,daysintotal*24,periods+1)]
        # Fills out the idx array, mapping each start index to an array of hours for that entire period
        hourlyarray = [collect(hourlyidx[i-1]+1:hourlyidx[i]) for i in eachindex(hourlyidx)[2:end]]
        return hourlyarray
    else
        # Finds the indices of the start of each of the periods for a daily resolution
        dailyidx = [round(Int64, i) for i in range(1,daysintotal+1,periods+1)]
        # Fills out the idx array, mapping each start index to an array of days for that entire period
        dailyarray = [collect(dailyidx[i-1]:dailyidx[i]-1) for i in eachindex(dailyidx)[2:end]]
        # Maps each day to an array of hours for that day
        if daily
            return dailyarray
        else 
            return [vcat([collect(24*(day-1)+1:24*(day-1)+24) for day in days]...) for days in dailyarray]
        end
    end
end

function normalizerows(
    matrix
)
    rowmatrix = repeat(sum(matrix, dims=1), size(matrix, 1), 1)
    matrix = matrix./rowmatrix
    return matrix
end

function orderingshuffle(
    ordering,
    std;
    it=100
)
    mean = 1
    lb, ub = 0, mean*2
    d = Truncated(Normal(mean, std*mean), lb, ub)
    shuffle = rand(d, size(ordering))
    for i in 1:it
        ordering = shuffle .* ordering
    end
    return normalizerows(ordering)
end

function extracttimeseriesdata(
    timeseries:: Timeseries_Data
)
    toi = (timeseries.solar, timeseries.wind_on, timeseries.wind_off, timeseries.price_DAM)
    return (
        mean(
            ts, 
            weights(repeat(decisionarray(timeseries.ordering), inner=24))) 
        for ts in toi
    )
end

function decisionarray(
    ordering::AbstractMatrix
)
    return vec(sum(ordering, dims=2))
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

- `basemodel::JuMP.Model`: Solved model with the results for the given parameters
"""
function solveoptimizationmodel(
    parameters::Parameter_Data, 
    time_limit::Int64,
    mip_gap::Float64
)
    # Defining the model
    basemodel = Model(Gurobi.Optimizer)

    # Setting the time limit and the mip gap
    set_time_limit_sec(basemodel, time_limit)
    set_optimizer_attribute(basemodel, "MIPgap", mip_gap)
    set_optimizer_attribute(basemodel, "DualReductions", 0)

    # Defining the most commonly used parameternames for readability
    scenario = parameters.scenario # Scenario data object
    hourly_target = scenario.hourly_target # Hourly target of hydrogen production in kg
    no_etc_periods = scenario.no_etc_periods # Number of periods for the electrical temporal correlation
    no_htc_periods = scenario.no_htc_periods # Number of periods for the hydrogen temporal correlation
    MPLDSP = scenario.MPLDSP # Minimum Partial Load of the DSP
    CFDSP = scenario.CFDSP # Capacity factor of the DSP
    EI_limit = scenario.EI_limit # Maximum EI in gCO2/kgH2
    lambda_TSO = scenario.lambda_TSO # Tariff for energy purchased from the grid

    components = parameters.components
    electrolyzer = components["electrolyzer"]
    storage = components["storage"]
    compressor = components["compressor"]
    
    timeseries = parameters.timeseries # Timeseries data object
    decision_array = repeat(decisionarray(timeseries.ordering), inner=24) # Array of weights for the representative days
    representative_days = size(timeseries.ordering, 1) # Number of representative days
    representative_timesteps = representative_days*24 # Number of timesteps in the representative days
    total_days = size(timeseries.ordering, 2) # Number of days in the total timeseries
    total_timesteps = total_days*24  # Number of timesteps in the total timeseries
    total_production = hourly_target*total_timesteps # Total hydrogen production in kg

    htc_isdaily = false#(no_htc_periods<=total_days) # Boolean for daily hydrogen temporal correlation, speed ups are iMPLDSPemented if true
    htc_periods = periodsarray(no_htc_periods, daysintotal=total_days, daily=htc_isdaily) # Array with the periods for the hydrogen temporal correlation
    etc_isdaily = false#(no_etc_periods<=total_days) # Boolean for daily electrical temporal correlation, speed ups are iMPLDSPemented if true

    if !isnothing(no_etc_periods) etc_periods = periodsarray(no_etc_periods, daysintotal=total_days, daily=etc_isdaily) end# Array with the periods for the electrical temporal correlation

    hourlymatching = (no_etc_periods==8760) # Boolean for hourly matching, speed ups are iMPLDSPemented if true

    # Defining the variables
        # Design (these caps are currently arbitrary, should be added to json file)
    @variables(basemodel, begin
        0 <= c_solar #Solar capacity
        0 <= c_wind_on #Onshore capacity
        0 <= c_wind_off #Offshore capacity
        0 <= c_electrolyzer #<= maxcapelectrolyzer(parameters) #Electrolyzer capacity
        0 <= c_compressor #Compressor capacity
        0 <= c_storage #<= maxcapstorage(parameters) #Storage capacity
    end)

    # Then each individual generation capacity, if given
    for (string, var) in zip(("solar", "wind_on", "wind_off"), (c_solar, c_wind_on, c_wind_off))
        if !isa(components[string].maxcapacity, Bool)
            @constraint(basemodel, var <= components[string].maxcapacity)
        elseif !isa(scenario.maxcapgeneration, Bool)
            @constraint(basemodel, var <= scenario.maxcapgeneration)
        end
    end

    # Optional constraints for fixed capacities
    #  @constraints(basemodel, begin
    #     c_storage == 10000
    # #     c_solar == 10
    # #     c_wind_on == 0
    # #     c_wind_off == 30
    #     c_electrolyzer == 50
    #end)

        # Operation
    @variables(basemodel, begin
        0 <= p_DAM_sell[1:representative_timesteps] #Power sold to the DAM in MW
        0 <= p_DAM_buy[1:representative_timesteps] # Power bought from the DAM in MW
        0 <= p_electrolyzer[1:representative_timesteps] #Electrolyzer power in MW
        0 <= p_compressor[1:representative_timesteps] #Compressor power in MW
        os_electrolyzer[1:representative_timesteps], Bin #Electrolyzer on/standby status
        0 <= h_electrolyzer[1:representative_timesteps] #Electrolyzer hydrogen production in kg
        0 <= h_storage_soc[1:total_timesteps] #Electrolyzer hydrogen storage state of charge in kg
        0 <= h_storage_in[1:representative_timesteps] #Electrolyzer hydrogen storage input in kg
        0 <= h_storage_out[1:representative_timesteps] #Electrolyzer hydrogen storage output in kg
        0 <= h_demand_dir[1:representative_timesteps] #Hydrogen output directly in kg
        0 <= h_demand[1:representative_timesteps] #Hydrogen output in kg
    end)

    # Defining the constraints

    # If hourly matching, we limit the potential power bought from the DAM s.t. we can't use it for electrolysis
    if hourlymatching
        @constraint(
            basemodel,
            hourlymatching_powercap_con[t in 1:representative_timesteps],
            p_DAM_buy[t] <= (1-os_electrolyzer[t])*c_electrolyzer*electrolyzer.P_standby+p_compressor[t]
        )
    end

    if !isnothing(no_etc_periods)
        # Defining the net consumption as the power bought, minus all the power sold and used for non-electrolyis purposes
        netconsumption_rep = [
            p_DAM_buy[t] - (
                p_DAM_sell[t] +
                p_compressor[t] +
                ((1 - os_electrolyzer[t]) * (electrolyzer.P_standby * c_electrolyzer)))
            for t in 1:representative_timesteps
        ]

        # Extending the net consumption to the total timeseries
        netconsumption_ext = rep_to_ext(netconsumption_rep, timeseries.ordering, daily=etc_isdaily)

            # Electrical temporal correlation constraint
        @constraint(
            basemodel, 
            etc_con[period in etc_periods], 
            (sum(
                netconsumption_ext[t]
            for t in period) <= 0)
        )
    end
    
        # Number of upper bounds on variables that increase the basemodel speed
            # Upper bound electrolyzer capacity, which should never exceed the max rated power of the RES
    @constraint(
        basemodel, 
        electrolyzer_c_upper_con,
        (c_electrolyzer <= (c_wind_on+c_wind_off+c_solar))
    )

        # Physical constraints
            # Energy balance
    @constraint(
        basemodel, 
        e_balance_con[t in 1:representative_timesteps], 
        (p_electrolyzer[t] + p_compressor[t] + p_DAM_sell[t] == (
            c_solar*timeseries.solar[t] +
            c_wind_on*timeseries.wind_on[t] +
            c_wind_off*timeseries.wind_off[t] + 
            p_DAM_buy[t]))
    )

        #  Upper bound electrolyzer power
    @constraint(
        basemodel, 
        electrolyzer_p_upper_con[t in 1:representative_timesteps], 
        (p_electrolyzer[t] <= (
        c_electrolyzer*os_electrolyzer[t] +
        electrolyzer.P_standby*c_electrolyzer*(1-os_electrolyzer[t])))
    )

        #  Lower bound electrolyzer power
    @constraint(
        basemodel, 
        electrolyzer_p_lower_con[t in 1:representative_timesteps], 
        (p_electrolyzer[t] >= (
        c_electrolyzer*electrolyzer.P_min*os_electrolyzer[t] +
        electrolyzer.P_standby*c_electrolyzer*(1-os_electrolyzer[t])))
    )

        # Electrolyzer hydrogen production
    @constraint(
        basemodel,
        electrolyzer_h_production_con[t in 1:representative_timesteps],
        (h_electrolyzer[t] == (
            electrolyzer.efficiency * 
            p_electrolyzer[t] * 
            os_electrolyzer[t]))
    )

    h_storage_in_ext = rep_to_ext([h_storage_in[t] for t in 1:representative_timesteps], timeseries.ordering)
    h_storage_out_ext = rep_to_ext([h_storage_out[t] for t in 1:representative_timesteps], timeseries.ordering)

        # Hydrogen storage capacity
    @constraint(
        basemodel,
        storage_capacity_con[t in 1:total_timesteps],
        (h_storage_soc[t] <= c_storage)
    )
    
        # Hydrogen storage initialization
    @constraint(
        basemodel,
        storage_soc_init_con,
        (h_storage_soc[1] == storage.soc_init*c_storage + h_storage_in_ext[1] - h_storage_out_ext[1])
    )
    
        # Electrolyzer hydrogen storage state of charge
    @constraint(
        basemodel,
        storage_soc_con[t in 2:total_timesteps],
        (h_storage_soc[t] == h_storage_soc[t-1] + h_storage_in_ext[t] - h_storage_out_ext[t])
    )

        # Hydrogen storage final state of charge equal to initial state of charge
    @constraint(
        basemodel,
        storage_conservation_con,
        h_storage_soc[total_timesteps] >= storage.soc_init*c_storage
    )

    # Compressor power consumption
    @constraint(
        basemodel,
        compressor_consumption_con[t in 1:representative_timesteps],
        p_compressor[t] == compressor.constant * h_storage_in[t]
    )

    @constraint(
        basemodel,
        compressor_capacity_con[t in 1:representative_timesteps],
        p_compressor[t] <= c_compressor
    )
   
    # Electrolyzer output
    @constraint(
        basemodel,
        electrolyzer_output_con[t in 1:representative_timesteps],
        (h_electrolyzer[t] == h_demand_dir[t] + h_storage_in[t])
    )

    # System output
    @constraint(
        basemodel,
        output_con[t in 1:representative_timesteps],
        (h_demand[t] == h_demand_dir[t] + h_storage_out[t])
    )

    #    Hourly flexibility constraint
    @constraint(
        basemodel,
        flex_lb_con[t in 1:representative_timesteps],
        MPLDSP*hourly_target/CFDSP <= h_demand[t]
    )

    @constraint(
        basemodel,
        flex_ub_con[t in 1:representative_timesteps],
        h_demand[t] <= hourly_target/CFDSP
    )

    #Hydrogen temporal correlation constraint
    h_demand_ext = rep_to_ext([h_demand[t] for t in 1:representative_timesteps], timeseries.ordering, daily=htc_isdaily)

    @constraint(
        basemodel, 
        htc_con[period in 1:length(htc_periods)],
        (sum(h_demand_ext[t] for t in htc_periods[period]) >= 
        length(htc_periods[period])*hourly_target)
    )

    if !isnothing(EI_limit)
        # Emission intensity constraint
        @constraint(
            basemodel,
            emission_con,
            (sum(decision_array[t]*p_DAM_buy[t]*timeseries.EI[t] for t in 1:representative_timesteps)/total_production <= 0.001*EI_limit)
        )
    end

    # Defining the objective function
    @objective(
        basemodel, 
        Min, 
        (
            getinvestmentcosts(basemodel, parameters) +
            sum(
                decision_array[t] * (
                p_DAM_buy[t] * (timeseries.price_DAM[t] + lambda_TSO) - 
                p_DAM_sell[t] * timeseries.price_DAM[t])
                for t in 1:representative_timesteps)
        )
    )

    # Solving the basemodel
    optimize!(basemodel)

    return basemodel
end

function getopportunity(
    parameters::Parameter_Data;
    time_limit::Int64=60
)

    # Defining some variables for visibility
    scenario = parameters.scenario
    components = parameters.components
    timeseries = parameters.timeseries

    decision_array = repeat(decisionarray(timeseries.ordering), inner=24)
    representative_timesteps = size(timeseries.ordering, 1)*24

    # Defining the model
    RESmodel = Model(Gurobi.Optimizer)

    # Setting the time limit and the mip gap
    set_time_limit_sec(RESmodel, time_limit)
    # set_optimizer_attribute(RESmodel, "DualReductions", 0)

    @variables(RESmodel, begin
        0 <= c_solar
        0 <= c_wind_on
        0 <= c_wind_off
    end) #Capacity of the RES

    capacityarray = [c_solar, c_wind_on, c_wind_off]
    componentarray = [components["solar"], components["wind_on"], components["wind_off"]]

    # Limiting individual capacities
    for (comp, var) in zip(componentarray, capacityarray)
        if !isa(comp.maxcapacity, Bool)
            @constraint(RESmodel, var <= comp.maxcapacity)
        elseif !isa(scenario.maxcapgeneration, Bool)
            @constraint(RESmodel, var <= scenario.maxcapgeneration)
        end
    end

    RESgen =  sum(cap * cf for (cap, cf) in zip(
                capacityarray, (timeseries.solar, timeseries.wind_on, timeseries.wind_off)))

    excludeinvest = ["electrolyzer", "storage", "compressor"]

    @objective(
        RESmodel,
        Min,
        (getinvestmentcosts(RESmodel, parameters, except=excludeinvest)- 
        sum(
            decision_array[t] * RESgen[t] * timeseries.price_DAM[t]
            for t in 1:representative_timesteps))
    )

    optimize!(RESmodel)

    return objective_value(RESmodel)

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
    verbose::Bool=false,
    opportunitycost::Number = 0
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
    total_days = size(parameters.timeseries.ordering, 2)
    total_timesteps = total_days*24
    decision_array = repeat(decisionarray(parameters.timeseries.ordering), inner=24)
    timeseries = parameters.timeseries
    ordering = timeseries.ordering
    price_DAM = timeseries.price_DAM

    # Defining the variables for the timeseries
    REScaps = (:c_solar, :c_wind_on, :c_wind_off)
    RESpdict = Dict(:c_solar=>:p_solar, :c_wind_on=>:p_wind_on, :c_wind_off=>:p_wind_off)
    REScf = (timeseries.solar, timeseries.wind_on, timeseries.wind_off)
    pchars = (:p_electrolyzer, :p_compressor, :p_DAM_buy, :p_DAM_sell)
    Hcaps = (:c_electrolyzer, :c_storage, :c_compressor)
    hchars = (:h_electrolyzer, :h_storage_in, :h_storage_out, :h_demand, :h_demand_dir)

    netconsumption = (
        var_data[:p_DAM_buy] .- (
            var_data[:p_DAM_sell] .+
            var_data[:p_compressor] .+
            ((1 .- var_data[:os_electrolyzer]) .* (parameters.components["electrolyzer"].P_standby * var_data[:c_electrolyzer])))
    )

    # Creating the timeseries dataframe
    timeseries_data = Dict()
    
    timeseries_data[:repdata] = DataFrame(
        merge(
            Dict(:datetime => datetimefromrepdays(timeseries.representativedaysarray)),
            # Dict with RES data
            Dict(RESpdict[REScap] => (ts*var_data[REScap]) for (ts, REScap) in zip(REScf, REScaps)),
            # Dict with other power data
            Dict(pchar => var_data[pchar] for pchar in pchars),
            # Dict with net consumption
            Dict(:netconsumption => netconsumption),
            # Dict with hydrogen data
            Dict(hchar => var_data[hchar] for hchar in hchars),
            # Dict with the day-ahead price
            Dict(:price_DAM => timeseries.price_DAM),
            # Dict with the emission data
            Dict(:EI => timeseries.EI, :attremissions => var_data[:p_DAM_buy].*timeseries.EI),
            # Dict with decision_array
            Dict(:decision_array => decision_array)
        )
    )
    timeseries_data[:extdata] = DataFrame(
        :datetime => datetimewholeyear(year, total_days),
        # Hydrogen storage state of charge
        :h_storage_soc => var_data[:h_storage_soc],
    )

    capacitycomponentpairs_RES = zip(
        (var_data[char] for char in REScaps),
        (parameters.components[component] for component in ("solar", "wind_on", "wind_off"))
    )
    capacitycomponentpairs_H = zip(
        (var_data[char] for char in Hcaps),
        (parameters.components[component] for component in ("electrolyzer", "storage", "compressor"))
    )

    total_cost = objective_value(model)-opportunitycost
    total_production = (parameters.scenario.hourly_target*total_timesteps)
    total_cost_DAM = sum(decision_array.*((var_data[:p_DAM_buy].*(price_DAM.+parameters.scenario.lambda_TSO)) .- var_data[:p_DAM_sell].*price_DAM))
    total_emissions = sum(decision_array.*var_data[:p_DAM_buy].*timeseries.EI)
    capital_cost_RES = sum(getcomponentcosts(capacity, parameters.scenario, component) for (capacity, component) in capacitycomponentpairs_RES)
    adjusted_cost_RES = (capital_cost_RES+total_cost_DAM-opportunitycost)
    RES_residual_factor = adjusted_cost_RES/capital_cost_RES
    LCOH = (total_cost)/total_production
    EI_H = total_emissions/total_production

    # Extracting outcomes
    FinancialDict = Dict(
        #Manual calculations
        "total_cost" => total_cost,
        "LCOH" => LCOH,
        "LCOH(unadjusted)" => objective_value(model)/(parameters.scenario.hourly_target*total_timesteps),
        "electrolyzer_cf" => mean((var_data[:p_electrolyzer]./var_data[:c_electrolyzer]), weights(decision_array)),
        "grid_cost" => total_cost_DAM,
        "grid_cost_fraction" => total_cost_DAM/total_cost,
        "grid_cost_LCOH" => total_cost_DAM/total_production,
        "total_emissions" => total_emissions,
        "EI_H" => EI_H,
        "AAC" => (LCOH-1.15)/(9.28-EI_H),
        "RES_capex_LCOH" => capital_cost_RES/total_production,
        "RES_cost" => adjusted_cost_RES,
        "RES_fraction" => adjusted_cost_RES/total_cost,
        "opportunity_cost" => opportunitycost,
    )

    println("opportunitycost: ", opportunitycost)
    
    RESDict = Dict()
    for (capacity, component) in capacitycomponentpairs_RES
        capex = getcomponentcosts(capacity, parameters.scenario, component)
        adjustedcost = capex*RES_residual_factor
        merge!(
            RESDict,
            Dict(
                "$(name(component))_capacity"=> capacity, 
                "$(name(component))_cost" => adjustedcost, 
                "$(name(component))_capex" => capex, 
                "$(name(component))_capex_LCOH" => capex/total_production, 
                "$(name(component))_fraction" => adjustedcost/total_cost,
                "$(name(component))_LCOH" => adjustedcost/total_production
            )
        )
    end

    HDict = Dict()
    for (capacity, component) in capacitycomponentpairs_H
        capex = getcomponentcosts(capacity, parameters.scenario, component)
        merge!(
            HDict,
            
            Dict(
                "$(name(component))_capacity"=> capacity, 
                "$(name(component))_cost" => capex, 
                "$(name(component))_fraction" => capex/total_cost,
                "$(name(component))_LCOH" => capex/total_production
            )
        )
    end

    HDict["storage_system_cost"] = HDict["storage_cost"]+HDict["compressor_cost"]
    HDict["storage_system_fraction"] = HDict["storage_system_cost"]/total_cost
    HDict["storage_system_LCOH"] = HDict["storage_system_cost"]/total_production

    println(HDict)

    outcome_data = merge(FinancialDict, RESDict, HDict)

    # Printing the results if verbose
    if verbose
        print(json(outcome_data,4))
    end

    return Dict(
        :outcome_data => outcome_data,
        :timeseries_data => timeseries_data
    )
end

function addoutcometodf(
    results::DataFrame;
    outcomes::Vector{Any}=[]
)
    if isempty(outcomes)
        outcomes = keys(results[4, :results][:outcome_data])
    end
    for outcome in outcomes
        println("Adding outcome: $(outcome)")
        data=[]
        for result in results[:, "results"]
            if result isa AbstractDict
                push!(data, result[:outcome_data][outcome])
            elseif  !(result isa Number) || isnan(result) 
                push!(data, NaN)
            end
        end
        dictdata = filter(t -> (t isa AbstractDict), data)
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
    title::String="LCOH for various scenarios",
    C::Union{Symbol,String}=:LCOH,
    Clabel::String="LCOH [€/kg]",
    normalized::Bool=false,
    unitmap::Float64=1.0
)
    X=:no_etc_periods
    Y2=:MPLDSP
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
    noY2labels = length(Y2labels)

    # Creating the value arrays
    Xvals = repeat(range(1, Xscenarios), inner=convert(Int64,Yscenarios))
    Yvals = repeat(range(1, Yscenarios), convert(Int64,Xscenarios))
    Cvals = [convert(Float32, x) for x in (results[:, C])]

    # Normalizing the values if required
    if unitmap != 1
        Cvals = Cvals.*unitmap
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
        title=title,
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
        ylabel="MPL (%)",
        yticks=(
            range(
                start=1/noY2labels,
                stop=10-1/noY2labels,
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
    results::DataFrame,
    savepath::String,
    parameters::Parameter_Data
    )
    writeable = select(results, Not("results"))
    CSV.write("$(savepath)/variabledata_$(Dates.format(now(), "mmddTHHMM")).csv", writeable)
end

function saveresults(
    results::Dict,
    savepath::String,
    parameters::Parameter_Data
    )
    timedict = Dict(
        1 => "Y",
        4 => "Q",
        12 => "M",
        26 => "bW",
        52 => "W",
        365 => "D",
        8760 => "H"
    )
    scenario = "$(timedict[parameters.scenario.no_etc_periods])$(timedict[parameters.scenario.no_htc_periods])$(round(Int64,(parameters.scenario.MPLDSP*100)))"
    open("$(savepath)/timeseriesdata_$(scenario)_$(Dates.format(now(), "mmddTHHMM")).csv", "w") do f
        JSON.print(f,results[:outcome_data]) 
    end    
    CSV.write("$(savepath)/timeseriesrepdata_$(scenario)_$(Dates.format(now(), "mmddTHHMM")).csv", results[:timeseries_data][:repdata])
    CSV.write("$(savepath)/timeseriesextdata_$(scenario)_$(Dates.format(now(), "mmddTHHMM")).csv", results[:timeseries_data][:extdata])
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
    opportunity::Bool=true
)
    if verbose
        println("Parameters: ", parameters)
    end
    solvedmodel = solveoptimizationmodel(parameters, 4000, 0.05)
    if opportunity
        opportunitycost = getopportunity(parameters)
        if opportunitycost>0
            throw(ArgumentError("Opportunity cost is positive: $(opportunitycost)"))
        end
    end
    if termination_status(solvedmodel) == MOI.OPTIMAL
        println("Optimal solution found.")
        results = getresults(
            solvedmodel, 
            parameters, 
            verbose=verbose, 
            opportunitycost=opportunitycost)
    else
        println("No optimal solution found: $(termination_status(solvedmodel)) \nReturning model.")
        return solvedmodel
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
    opportunity::Bool=true
)
    if verbose
        println("Parameters: ", parameterdf)
    end
    solutiondataframe = select(parameterdf, Not("parameters"))
    results = []
    for (i, parameters) in enumerate(parameterdf[:, "parameters"])
        println("
        At $i / $(length(parameterdf[:,"parameters"]))
        ETC periods: $(parameters.scenario.no_etc_periods)
        HTC periods: $(parameters.scenario.no_htc_periods) 
        MPL: $(parameters.scenario.MPLDSP)")
        solvedmodel = solveoptimizationmodel(parameters, 4000, 0.05)
        if opportunity
            opportunitycost = getopportunity(parameters)
            if opportunitycost>0
                throw(ArgumentError("Opportunity cost is positive: $(opportunitycost)"))
            end
        end
        push!(
            results,
            getresults(
                solvedmodel, 
                parameters, 
                verbose=verbose,
                opportunitycost=opportunitycost)
        )
    end
    solutiondataframe[:, "results"] = results
    return solutiondataframe
end

function shufflescript(
    parameters::Parameter_Data,
    shuffles::Int64;
    verbose::Bool=true,
    std=0.25
)

    # Initiating dataframe to hold the results
    outcomedf = DataFrame(
        :CFsolar => [],
        :CFwind_on => [],
        :CFwind_off => [],
        :DAM => [],
        :Ordering => []
    )
    baseordering = copy(parameters.timeseries.ordering)
    for i in 1:shuffles
        println("Shuffle $i")
        # Shuffling the timeseries
        parameters = @set parameters.timeseries.ordering = orderingshuffle(baseordering, std)
        # Optimizing the model
        #result = optimizeptx(parameters, verbose=verbose, savecsv=false)

        delta = sum(abs.(parameters.timeseries.ordering-baseordering))
        # Adding the results to the dataframe
        push!(outcomedf, (extracttimeseriesdata(parameters.timeseries)..., parameters.timeseries.ordering))
    end
    return outcomedf
end


function mainscript(
    parameterfilename::String,
    verbose::Bool=true;
    savepath::String="$(home_dir)/Results",
    savecsv::Bool=false,
    savejld::Bool=true,
    opportunity::Bool=true
)   
    parameters = fetchparameterdata(parameterfilename)
    results = optimizeptx(parameters, verbose=verbose, savepath=savepath, opportunity=opportunity)
    if savecsv
        saveresults(results, savepath, parameters)
    end
    if savejld
        save_object("$(savepath)/results_$(Dates.format(now(), "mmddTHHMM")).jld2", results)
    end
    return results
end

