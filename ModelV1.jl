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
 

function optimizationmodel(parameters::PtxOpt.Parameter_Data)
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
        PtxOpt.lowerbound_electrolyzer(model, electrolyzer, k, t))
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
        (sum(PtxOpt.getcomponentcosts(parameters, (c_solar, c_wind, c_electrolyzer)) +
        sum(timeseries.price_DAM[k,t]*(p_DAM_buy[k,t]-p_DAM_sell[k,t]) for k in 1:tc_periods, t in 1:tc_length)))
    )

    # Solving the model
    optimize!(model)

    return model
end

# Printing the results
function showresults(model)
    results = getresultsdict(model)

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

    plot(range(1,simulation_length), [p_solar p_wind p_buy p_sell p_el], label=["Solar" "Wind" "Buy" "Sell" "Electrolyzer"], xlabel="Time", ylabel="Power (MW)", title="Power flow", legend=:topleft, legend_font_pointsize=legendfontsize)
    plot!(twinx(), price_subset, color=:black, linestyle=:dash, label="Price", ylabel="Price (€/MWh)", legend=:topright, show=true, legend_font_pointsize=legendfontsize)
    savefig("powerflow.svg")
end

"""
    main(parameterfilename::String)

"""
function main(parameterfilename::String)
    parameters = PtxOpt.fetchparameterdata(parameterfilename)
    model = optimizationmodel(parameters)
    if termination_status(model) == MOI.OPTIMAL
        println("Optimal solution found.")
        showresults(model)
    else
        throw(ErrorException("Optimal solution not found."))
    end
    return model
end