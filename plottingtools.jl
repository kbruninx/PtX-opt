using CSV
using DataFrames
using Statistics
using StatsBase
using CairoMakie
using Tables
using LinearAlgebra
using Dates
using KissSmoothing
using FileIO, JLD2
CairoMakie.activate!(type = "png")

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

# Loading data

# NLdata = ["NL60_resulting_profiles.csv", "NL60_ordering_variable.csv"]
# ESdata = ["ES60_resulting_profiles.csv", "ES60_ordering_variable.csv"]

NLMW80, NLHW80 = load_object("Results/NLMW80.jld2"), load_object("Results/NLHW80.jld2")
tsves = [NLHW80[:timeseries_data], NLMW80[:timeseries_data]]
# NLMW80GEO, NLHW80GEO = load_object("Results/NLMW80GEO.jld2"), load_object("Results/NLHW80GEO.jld2")
# tsgeo = [NLHW80GEO[:timeseries_data], NLMW80GEO[:timeseries_data]]

geostoragecompdf = load_object("Results/geostoragecompdf.jld2")
ves_ts = [row.results[:timeseries_data] for row in eachrow(sort(filter(row->row.label==:ves, geostoragecompdf), :no_etc_periods, rev=true))]
geo_ts = [row.results[:timeseries_data] for row in eachrow(sort(filter(row->row.label==:geo, geostoragecompdf), :no_etc_periods, rev=true))]

ves_out = [row.results[:outcome_data] for row in eachrow(sort(filter(row->row.label==:ves, geostoragecompdf), :no_etc_periods, rev=true))]
geo_out = [row.results[:outcome_data] for row in eachrow(sort(filter(row->row.label==:geo, geostoragecompdf), :no_etc_periods, rev=true))]

geobarplotdata = load_object("Results/geobarplot_data.jld2")

# orderingNL60 = transpose(CSV.File("NL60_ordering_variable.csv") |> Tables.matrix)

scenariosNL = load_object("Results/ScenariosNL(-HM25)full.jld2")

eidf = load_object("Results/EIcomparisondf.jld2")

# comparisonvarsMPL = [(8760,1,0.0), (8760,1,0.25), (8760,1,0.5), (8760,1,0.75), (8760,1,0.9)]
# comparisonscensMPL = filter(scenario -> (scenario[:no_etc_periods], scenario[:no_htc_periods], scenario[:MPLDSP]) in comparisonvarsMPL, scenariosNL)
# comparissontsMPL = [scenario[:Results][:timeseries_data] for scenario in eachrow(sort(comparisonscensMPL, [:MPLDSP]))]

# comparisonvarsHTC = [(8760,1,0.0), (8760,12,0.0), (8760,52,0.0), (8760,365,0.0), (8760,8760,0.0)]
# comparisonscensHTC = filter(scenario -> (scenario[:no_etc_periods], scenario[:no_htc_periods], scenario[:MPLDSP]) in comparisonvarsHTC, scenariosNL)
# comparissontsHTC = [scenario[:Results][:timeseries_data] for scenario in eachrow(sort(comparisonscensHTC, [:no_htc_periods]))]

# display(comparisonscensMPL)
# display(comparisonscensHTC)

# Auxiliary functions

function sum_periodically(
    profiledf::DataFrame;
    periods::Int64=52,
)
    daysintotal = Int64(nrow(profiledf)/24)
    # Finds the indices of the start of each of the periods for a daily resolution
    dailyidx = [round(Int64, i) for i in range(1,daysintotal+1,periods+1)]
    # Fills out the idx array, mapping each start index to an array of days for that entire period
    dailyarray = vcat([repeat([i], (dailyidx[i+1]-dailyidx[i])) for i in eachindex(dailyidx)[1:end-1]]...)
    hourlyarray = repeat(dailyarray, inner=24)
    profiledf[:, :period] = hourlyarray

    # Summing the data by the period
    periodicdf = select(combine(groupby(profiledf, :period), names(profiledf) .=> mean, renamecols=false), Not(:period))

    return periodicdf
end

function sum_periodically(
    profilevector::AbstractVector;
    periods::Int64=52,
)
    daysintotal = Int64(length(profilevector)/24)
    # Finds the indices of the start of each of the periods for a daily resolution
    dailyidx = [(24*round(Int64, i))+1 for i in range(0,daysintotal,periods+1)]
    # Fills out the idx array, mapping each start index to an array of days for that entire period
    periodicvec = [mean(profilevector[dailyidx[i]:dailyidx[i+1]-1]) for i in eachindex(dailyidx)[1:end-1]]
    return periodicvec
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



# PLOTS FOR THE REPORT IN ORDER OF APPEARANCE

# REPRESENTATIVE DAYS SUMMER & WINTER

function plotoperations(
    timeseriesdata::Dict,
)
    repdata = timeseriesdata[:repdata]
    extdata = timeseriesdata[:extdata]

    summerdays = [164,165]
    winterdays = [345,346]

    summerdata = filter(row -> dayofyear(row.datetime) in summerdays, repdata)
    summerdata = innerjoin(summerdata, extdata, on=:datetime)
    winterdata = filter(row -> dayofyear(row.datetime) in winterdays, repdata)
    winterdata = innerjoin(winterdata, extdata, on=:datetime)

    display(summerdata)
    display(winterdata)

    fig = Figure(resolution=(1200,700))
    
    # Creating the grid for the plots
    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()
    gc = fig[2,1] = GridLayout()
    gd = fig[2,2] = GridLayout()

    # Creating the axes
    axamain = Axis(ga[1,1], xgridvisible=false, ygridvisible=false)
    axadual = Axis(ga[1,1], yaxisposition=:right, xgridvisible=false, ygridvisible=false)
    axbmain = Axis(gb[1,1], xgridvisible=false, ygridvisible=false)
    axbdual = Axis(gb[1,1], yaxisposition=:right, xgridvisible=false, ygridvisible=false)
    axcmain = Axis(gc[1,1], xgridvisible=false, ygridvisible=false)
    axcdual = Axis(gc[1,1], yaxisposition=:right, xgridvisible=false, ygridvisible=false)
    axdmain = Axis(gd[1,1], xgridvisible=false, ygridvisible=false)
    axddual = Axis(gd[1,1], yaxisposition=:right, xgridvisible=false, ygridvisible=false)


    # Linking the axes of the same plots
    #linkxaxes!(axamain, axadual, axbmain, axbdual, axcmain, axcdual, axdmain, axddual)
    linkxaxes!(axamain, axadual, axcmain, axcdual)
    linkxaxes!(axbmain, axbdual, axdmain, axddual)
    linkyaxes!(axamain, axbmain)
    linkyaxes!(axadual, axbdual)
    linkyaxes!(axcmain, axdmain)
    linkyaxes!(axcdual, axddual)

    # hiding the x-axes of the top plots
    for ax in [axamain, axbmain]
        hidexdecorations!(ax, ticks=false)
    end

    # hiding the inside y-axes
    for ax in [axadual, axbmain, axcdual, axdmain]
        hideydecorations!(ax, ticks=false)
    end

    # hiding the dual x-axes
    for ax in [axadual, axbdual, axcdual, axddual]
        hidexdecorations!(ax)
    end

    # Titles for the rows
    axamain.title = "Example Period in Summer"
    axbmain.title = "Example Period in Winter"

    # Labeling the axes
    axamain.ylabel = "Power [MW]"
    axbdual.ylabel = "DAM price [€/MWh]"
    axcmain.ylabel = "Hydrogen massflow [kg/h]"
    axddual.ylabel = "Storage SOC [kg]"

    # Row 1 - Power flows

    labels1 = ["RES", "Electrolyzer", "Grid (bought)"]
    colors1 = [:gold, :dodgerblue3, :firebrick]

    # labels1 = ["RES", "Electrolyzer", "Grid (bought)", "Grid (sold)"]
    # colors1 = [:gold, :dodgerblue3, :mediumseagreen, :firebrick]
    datapoints1 = ["p_RES", "p_electrolyzer", "p_DAM_buy"]
    colortwin = :black

    # Grid A - Power flows summerdata

    summerdata[:, "p_RES"] = summerdata[:, "p_solar"] + summerdata[:, "p_wind_on"] + summerdata[:, "p_wind_off"]
    # summerpowerflows = [summerdata[:, datapoint] for datapoint in ["p_RES", "p_electrolyzer", "p_DAM_buy"]]

    summerpowerflows = [summerdata[:, datapoint] for datapoint in datapoints1]

    series!(axamain, summerpowerflows, labels=labels1, color=colors1)

    lines!(axadual, summerdata[:, "price_DAM"], label="DAM price", color=colortwin, linestyle=:dash)

    # Grid B - Power flows winterdata

    winterdata[:, "p_RES"] = winterdata[:, "p_solar"] + winterdata[:, "p_wind_on"] + winterdata[:, "p_wind_off"]
    # winterpowerflows = [winterdata[:, datapoint] for datapoint in ["p_RES", "p_electrolyzer", "p_DAM_buy"]]

    winterpowerflows = [winterdata[:, datapoint] for datapoint in datapoints1]

    series!(axbmain, winterpowerflows, labels=labels1, color=colors1)

    lines!(axbdual, winterdata[:, "price_DAM"], label="DAM price", color=colortwin, linestyle=:dash)

    # Legend 

    Legend(fig[1,3], axamain, valign=:top)
    Legend(fig[1,3], axadual, valign=:bottom)

    # Row 2 - Hydrogen flows

    labels2 = ["Electrolyzer", "DSP feed"]
    colors2 = [:dodgerblue3, :mediumseagreen]
    datapoints2 = ["h_electrolyzer", "h_demand"]
    colortwin = :black

    # Grid C - Hydrogen flows summerdata

    summerhydrogenflows = [summerdata[:, datapoint] for datapoint in datapoints2]
    
    series!(axcmain, summerhydrogenflows, labels=labels2, color=colors2)

    lines!(axcdual, summerdata[:, "h_storage_soc"], label="Storage SOC", color=colortwin, linestyle=:dash)

    # Grid D - Hydrogen flows winterdata

    winterhydrogenflows = [winterdata[:, datapoint] for datapoint in datapoints2]
    
    series!(axdmain, winterhydrogenflows, labels=labels2, color=colors2)

    lines!(axddual, winterdata[:, "h_storage_soc"], label="Storage SOC", color=colortwin, linestyle=:dash)

    # Formatting
    axamain.xticks = (range(1,48,5), repeat([""], 5))
    axbmain.xticks = (range(1,48,5), repeat([""], 5))
    axcmain.xticks = (range(1,48,5), ["Jun 12\n2020", "", "Jun 13", "", "Jun 14"])
    axdmain.xticks = (range(1,48,5), ["Dec 10\n2020", "", "Dec 11", "", "Dec 12"])

    round_step(x, step) = round(x / step) * step

    xlims!.([axadual, axbdual, axcdual, axddual], 1, 48)
    ylims!.([axamain, axbmain], low=-.7)
    ylims!.([axadual, axbdual], low=-5)

    powerflowticks = range(extrema(vcat(summerpowerflows..., winterpowerflows...))..., 5)
    dampriceticks = range(extrema(vcat(summerdata[:, "price_DAM"], winterdata[:, "price_DAM"]))..., 5)

    hydrogenflowticks = range(extrema(vcat(summerhydrogenflows..., winterhydrogenflows...))..., 5)

    # Legend

    Legend(fig[2,3], axcmain, valign=:top)
    Legend(fig[2,3], axcdual, valign=:bottom)

    for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end


    fig
end

# HEATMAP AND CONTOUR PLOTS

"""
The plotheatmap() function takes a DataFrame with the results of a scenario 
analysis and plots a heatmap of the results. The DataFrame should have the 
following columns:

- no_etc_periods: The number of periods in the ETC
- no_htc_periods: The number of periods in the HTC
- MPLDSP: The MPL of the DSP

The function also takes the following keyword arguments:

- title: The title of the plot
- C: The column to plot
- Clabel: The label of the colorbar
- normalized: Whether to normalize the values
- unitmap: A factor to multiply the values with
"""
function plotheatmap(
    results::DataFrame,
    gl::GridLayout;
    C::Union{Symbol,String}=:LCOH,
    colorrange=nothing,
    normalized::Bool=false,
    unitmap::Float64=1.0
)
    X=:no_etc_periods
    Y2=:MPLDSP
    Y1=:no_htc_periods


    # Sorting results
    sort!(results, [X, Y1, Y2])

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
        Cvals = (Cvals.-mean(filter(!isnan, Cvals)))./mean(filter(!isnan, Cvals)).*100
        println(Cvals)
    end

    ax1=Axis(
        gl[1,1], 
        xlabel="ETC length", 
        ylabel="HTC length", 
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

    if isnothing(colorrange)
        hm = heatmap!(Xvals, Yvals, Cvals, colormap=cgrad(:RdYlGn_8, rev=true), nan_color = :snow2)
    else
        hm = heatmap!(Xvals, Yvals, Cvals, colormap=cgrad(:RdYlGn_8, rev=true), nan_color = :snow2, colorrange=colorrange)
    end

    ax2=Axis(
        gl[1,1],
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

    return ax1, ax2, hm
end

function plotcontour(
    results::DataFrame,
    gl::GridLayout;
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

    println(Xvals)
    println(Yvals)
    println(Cvals)

    # Normalizing the values if required
    if unitmap != 1
        Cvals = Cvals.*unitmap
    end
    if normalized
        Cvals = Cvals./mean(filter(!isnan, Cvals)).*100
        Clabel = "Deviation from the mean [%]"
    end

    ax1=Axis(
        gl[1,1], 
        xlabel="ETC length", 
        ylabel="HTC length", 
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
        yticksize=0,
    )

    #xlims!(ax1, 0.5, Xscenarios+0.5)

    hm = contourf!(Xvals, Yvals, Cvals, colormap=cgrad(:RdYlGn_8, rev=true), nan_color = :snow2)

    ax2=Axis(
        gl[1,1],
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

    tightlimits!(ax1)

    return ax1, ax2, hm
end

function plotLCOHhm(
    results
)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()

    ax1, ax2, hm = plotheatmap(results, ga, C=:LCOH)

    Colorbar(gb[1,1], hm, label="LCOH [€/kg]")

    ax1.title = "LCOH"

    fig
end


function plotEMhm(
    results
)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()

    ax1, ax2, hm = plotheatmap(results, ga, C=:average_emissions)

    Colorbar(gb[1,1], hm, label=rich("Attributable Emissions Grid Electricity [kg", subscript("CO2") ,"/kg",  subscript("H2"),")"))

    ax1.title = "Attributable emissions for various scenarios"

    fig
end

function plotAAChm(
    results
)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()

    ax1, ax2, hm = plotheatmap(results, ga, C=:AAC)

    Colorbar(gb[1,1], hm, label=rich("Attributable Abatement Costs [€/kg",  subscript("CO2"),")"))

    ax1.title = "Attributable abatement costs for various scenarios"

    fig
end

function plotLCOHcont(
    results
)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()

    ax1, ax2, hm = plotcontour(results, ga, C=:LCOH)

    Colorbar(gb[1,1], hm, label="LCOH [€/kg]")
    fig
end

function plotheatmapcomponents(
    results::DataFrame
)
    colorrange = (-200, 200)

    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()
    gc = fig[2,1] = GridLayout()
    gd = fig[2,2] = GridLayout()
    ge = fig[:,3] = GridLayout()

    axamain, axadual, hma = plotheatmap(results, ga, C=:RES_capital_cost, colorrange=colorrange, normalized=true)
    axbmain, axbdual, hmb = plotheatmap(results, gb, C=:electrolyzer_cost, colorrange=colorrange, normalized=true)
    axcmain, axcdual, hmc = plotheatmap(results, gc, C=:storage_system_cost, colorrange=colorrange, normalized=true)
    axdmain, axddual, hmd = plotheatmap(results, gd, C=:Grid_revenue, colorrange=colorrange, normalized=true)

    axamain.title = "RES"
    axbmain.title = "Electrolyzer"
    axcmain.title = "Storage system"
    axdmain.title = "Grid costs/profits"

    for ax in [axamain, axbmain]
        hidexdecorations!(ax, minorticks=false)
    end
    for ax in [axadual, axcdual, axbmain, axdmain]
        hideydecorations!(ax, ticks=false, minorticks=false)
    end

    Colorbar(ge[1,1], colormap=cgrad(:RdYlGn_8, rev=true), label="Fraction of total cost [%]", colorrange=colorrange)
    fig
end

# PRODUCTION LEVEL INTERMITTENCY FOR HTC AND MPL

function compareproductionlevels(
    timeseriesdfsMPL::Vector{Dict{Any, Any}},
    timeseriesdfsHTC::Vector{Dict{Any, Any}};
    labelsMPL::Vector{String}=["HY00", "HY25", "HY50", "HY75", "HY90"],
    labelsHTC::Vector{String}=["HY00", "HM00", "HW00", "HD00", "HH00"],
    ordering::AbstractMatrix=orderingNL60
)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()

    axa = Axis(ga[1,1], title="Weekly average hydrogen production for MPL scenarios", ylabel="Mean production per week [kg/h]", ylabelfont=:bold, xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xticksvisible=false)
    axb = Axis(gb[1,1], title="Weekly average hydrogen production for HTC scenarios", xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xticksvisible=false)
    colors = [:mediumseagreen, :dodgerblue3, :gold, :darkorange, :firebrick]

    msize = 5

    notimesteps = 52
    stepspermonth = notimesteps/12
    Xticks = collect(range(stepspermonth/2, notimesteps-(stepspermonth/2), 12))
    Xticklabels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for (i, timeseriesdf) in enumerate(timeseriesdfsMPL)
        hDSP = rep_to_ext(timeseriesdf[:repdata][:, :h_demand], ordering)
        hDSPweekly = sum_periodically(hDSP, periods=52)

        scatterlines!(axa, hDSPweekly, color=colors[i], label=labelsMPL[i], markersize=msize)

    end
    for (i, timeseriesdf) in enumerate(timeseriesdfsHTC)
        hDSP = rep_to_ext(timeseriesdf[:repdata][:, :h_demand], ordering)
        hDSPweekly = sum_periodically(hDSP, periods=52)

        scatterlines!(axb, hDSPweekly, color=colors[i], label=labelsHTC[i], markersize=msize)

    end

    Legend(fig[2,1], axa, orientation = :horizontal, tellwidth = false, tellheight = true)
    Legend(fig[2,2], axb, orientation = :horizontal, tellwidth = false, tellheight = true)

    axa.xticks = (Xticks, Xticklabels)
    axb.xticks = (Xticks, Xticklabels)

    xlims!(axa, 1, notimesteps)
    xlims!(axb, 1, notimesteps)

    for (label, layout) in zip(["A", "B"], [ga, gb])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    return fig
end

# BAR PLOTS

scenariosNL[:, "solar_capex_LCOH"] = [getcomponentcosts(row.solar_capacity, paramsNL.scenario, paramsNL.components["solar"])/(8760*100) for row in eachrow(scenariosNL)]
scenariosNL[:, "wind_on_capex_LCOH"] = [getcomponentcosts(row.wind_on_capacity, paramsNL.scenario, paramsNL.components["wind_on"])/(8760*100) for row in eachrow(scenariosNL)]

function makebarplotscenarios(
    results::DataFrame;
    includeres=false
)
    #Extracting the data
    components = ["grid_cost_LCOH", "wind_on_capex_LCOH", "solar_capex_LCOH", "electrolyzer_LCOH", "storage_system_LCOH", "LCOH"]
    variablesofinterest = ["no_etc_periods", "no_htc_periods", "MPLDSP"]
    voidict = Dict(
        "no_etc_periods" => "ETC length",
        "MPLDSP" => "MPL",
        "no_htc_periods" => "HTC length"
    )
    # Selecting the components of interest
    parsedresults = results[:, filter(name -> name in components || name in variablesofinterest, names(results))]

    # Preparing the figures
    fig = Figure(resolution=(1200,700))

    ga = fig[1:2,1] = GridLayout()
    gb = fig[1,2] = GridLayout()
    gc = fig[2,2] = GridLayout()

    # for comp in [:grid_cost, :RES_capital, :electrolyzer, :storage_system]
    #     sort!(parsedresults, comp)
    # end

    # constructing the axes
    axes1 = []
    # axes2 = []
    for (voi, pos) in zip(variablesofinterest, [ga[1,1], gb[1,1], gc[1,1]])
        push!(axes1, Axis(
                        pos,
                        xlabel=voidict[voi], 
                        ylabel="LCOH [€/kg]",
                        xminorticksvisible=true,
                        xticksvisible=false,
                        xgridvisible=false,
        ))
        # if includeres
        #     push!(axes2, Axis(
        #                     pos,
        #                     yaxisposition=:right,
        #                     ylabel="RES capacity [MW]",
        #                     xminorticksvisible=false,
        #                     xticksvisible=false,
        #                     xgridvisible=false,
        #                     ygridvisible=false,
        #     ))
        #     hidexdecorations!(axes2[end])
        #     hideydecorations!(axes2[end])
        # else
        #     push!(axes2, nothing)
        # end
    end

    for (i, voi) in enumerate(variablesofinterest)
        # Selecting the results for the variable of interest, getting rid of the other variables of interest
        nonvoi = filter(var -> var != voi, variablesofinterest)
        voidf = select(parsedresults, Not(nonvoi))
        # Grouping the results by the variable of interest
        
        meandf = combine(groupby(voidf, voi), Not(voi) .=> mean, renamecols=false)

        display(meandf)

        sort!(meandf, voi)

        xticklabels = unique(meandf[:, voi])
        if voi in ["no_etc_periods", "no_htc_periods"]
            xticklabels = [timedict[i] for i in xticklabels]
        else 
            xticklabels = ["$(100*label)%" for label in xticklabels]
        end
        
        meandict = Dict(
            "LCOH" => meandf[:, :LCOH],
            "grid_cost_LCOH" => meandf[:, :grid_cost_LCOH],
            "wind_on_capex_LCOH" => meandf[:, :wind_on_capex_LCOH],
            "solar_capex_LCOH" => meandf[:, :solar_capex_LCOH],
            "electrolyzer_LCOH" => meandf[:, :electrolyzer_LCOH],
            "storage_system_LCOH" => meandf[:, :storage_system_LCOH]
        )
        show(meandict)

        if includeres
            resdict = Dict(
                "solar_capacity" => meandf[:, :solar_capacity],
                "wind_on_capacity" => meandf[:, :wind_on_capacity]
            )
        else
            resdict = nothing
        end

        fig, axes1[i] = makebarplot(meandict, fig, axes1[i], xticklabels, resdict=resdict)
    end
    # final formatting 
    Legend(fig[3,:], axes1[1], orientation = :horizontal, tellwidth = false, tellheight = true)
    linkyaxes!(axes1...)
    # if includeres
    #     #Legend(fig[3,2], axes2[1], orientation = :horizontal, tellwidth = false, tellheight = true)
    #     linkyaxes!(axes2...)
    #     for i in eachindex(axes1)
    #         linkyaxes!(axes1[i], axes2[i])
    #     end
    # end

    for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    fig
end

labelEIcomparison = [
    rich("HW80\n(EI=17g", subscript("CO2") ,"/kg",  subscript("H2"),")"), 
    rich("EIW80\n(EI<=17g", subscript("CO2") ,"/kg",  subscript("H2"),")"),
    rich("MW80\n(EI=5kg", subscript("CO2") ,"/kg",  subscript("H2"),")"),
    rich("EIW80\n(EI=5kg", subscript("CO2") ,"/kg",  subscript("H2"),")")]

function makebarplotmanual(
    resultdicts::AbstractVector;
    labels=["HW80\nVessel", "HW80\nGeological", "MW80\nVessel", "MW80\nGeological"],
    includeres=false

)

    fig = Figure(resolution=(1400,700))

    ga = fig[1,1] = GridLayout()

    ax = Axis(
        ga[1,1],
        xlabel="Scenario", 
        xlabelfont=:bold,
        ylabel="LCOH [€/kg]",
        ylabelfont=:bold,
        xminorticksvisible=true,
        xticksvisible=false,
        xgridvisible=false,
        title="Bar plot of the LCOH contributions ETC and EI constraints"
    )

    symbols = ["LCOH", "grid_cost_LCOH", "wind_on_capex_LCOH", "solar_capex_LCOH", "electrolyzer_LCOH", "storage_system_LCOH"]

    jointdict = Dict(
        symb => vcat([resultdict[symb] for resultdict in resultdicts]...) for symb in symbols
    )

    if includeres
        resdict = Dict(
            "solar_capacity" => vcat([resultdict["solar_capacity"] for resultdict in resultdicts]...),
            "wind_on_capacity" => vcat([resultdict["wind_on_capacity"] for resultdict in resultdicts]...)
        )
    else
        resdict = nothing
    end

    fig, ax = makebarplot(jointdict, fig, ax, labels)

    Legend(fig[2,1], ax, orientation = :horizontal, tellwidth = false, tellheight = true)

    return fig
end

function makebarplot(
    dict::Dict,
    fig,
    ax,
    xticklabels;
    resdict=nothing
)
    includeres = !isnothing(resdict)
    noscen = length(xticklabels) # Number of scenarios
    #nocomp = length(components) # Number of components


    X = collect(1:noscen) # All the data (length nocomp X nodata)


    storage = dict["LCOH"]
    electrolyzer = dict["LCOH"].- dict["storage_system_LCOH"]
    solar = dict["solar_capex_LCOH"] .+ dict["wind_on_capex_LCOH"] .+ dict["grid_cost_LCOH"]
    wind_on = dict["wind_on_capex_LCOH"] .+ dict["grid_cost_LCOH"]
    grid = dict["grid_cost_LCOH"]

    Y = [storage, electrolyzer, solar, wind_on, grid]
    complabels = ["Storage", "Electrolyzer", "Solar", "Wind (onshore)", "Grid"]
    colors = [:firebrick, :gold, :lightskyblue1, :dodgerblue4, :mediumseagreen]

    display(Y)
    for i in eachindex(Y)
        barplot!(ax, X, Y[i], 
        color=colors[i], 
        label=complabels[i], 
        gap=0.1, 
        strokecolor=:black, 
        strokewidth=1,)
    end

    ax.xticks = (range(1,noscen), xticklabels)
    ax.xminorticks = range(0.5, stop=noscen+0.5, length=noscen+1)

    xlims!(ax, 0.5, noscen+0.5)

    return fig, ax
end


# Geocomparison chapter 4
function plotstoragelevels(tscomp, tsgeo, tslabels)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()
    gc = fig[2,1] = GridLayout()
    gd = fig[2,2] = GridLayout()

    axa = Axis(ga[1,1], title="Pressurized Vessel Storage", ylabel="Storage level [tonnes]", xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xticksvisible=false)
    axb = Axis(gb[1,1], title="Geological Storage", xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xticksvisible=false)
    axc = Axis(gc[1,1], ylabel="Monthly average RES production [MW]", xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xticksvisible=false)
    axd = Axis(gd[1,1], xgridvisible=false, ygridvisible=false, xminorticksvisible=true, xticksvisible=false)
    # axcdual = Axis(gc[1,1], yaxisposition=:right, yticksvisible=false, ygridvisible=false, yminorticksvisible=true)
    # axddual = Axis(gd[1,1], ylabel="Monthly average revenue [€/h]", yaxisposition=:right, ygridvisible=false)

    # hiding the x-axes of the top plots
    for ax in [axa, axb]
        hidexdecorations!(ax, minorticks=false)
    end

    ordering = transpose(CSV.File("NL60_ordering_variable.csv") |> Tables.matrix)

    colors = [:dodgerblue3, :gold, :firebrick, :mediumseagreen]
    selectedcolors = colors[1:length(tscomp)]

    # Collecting data
        # Storage
    storagesoc_comp = [ts[:extdata][:, :h_storage_soc]./1000 for ts in tscomp]
    storagesoc_comp_daily = sum_periodically.(storagesoc_comp, periods=365)
    storagesoc_comp_monthly = sum_periodically.(storagesoc_comp, periods=12)

    storagesoc_geo = [ts[:extdata][:, :h_storage_soc]./1000 for ts in tsgeo]
    storagesoc_geo_daily = sum_periodically.(storagesoc_geo, periods=365)
    storagesoc_geo_monthly = sum_periodically.(storagesoc_geo, periods=12)

        # RES
    res_comp = [rep_to_ext((ts[:repdata][:, :p_solar]+ts[:repdata][:, :p_wind_on]), ordering) for ts in tscomp]
    res_comp_monthly = sum_periodically.(res_comp, periods=12)
    res_comp_weekly = sum_periodically.(res_comp, periods=52)
    res_comp_daily = sum_periodically.(res_comp, periods=365)

    res_geo = [rep_to_ext((ts[:repdata][:, :p_solar]+ts[:repdata][:, :p_wind_on]), ordering) for ts in tsgeo]
    res_geo_monthly = sum_periodically.(res_geo, periods=12)
    res_geo_weekly = sum_periodically.(res_geo, periods=52)
    res_geo_daily = sum_periodically.(res_geo, periods=365)

        # Revenue DAM
    revenue_comp = [rep_to_ext(ts[:repdata][:, "p_DAM_sell"].*ts[:repdata][:, "price_DAM"], ordering) for ts in tscomp]
    revenue_comp_monthly = sum_periodically.(revenue_comp, periods=52)

    revenue_geo = [rep_to_ext(ts[:repdata][:, "p_DAM_sell"].*ts[:repdata][:, "price_DAM"], ordering) for ts in tsgeo]
    revenue_geo_monthly = sum_periodically.(revenue_geo, periods=52)


    for arrays in [storagesoc_comp_monthly, storagesoc_geo_monthly, res_comp_monthly, res_geo_monthly]
        push!(arrays[1], arrays[1][end])
        push!(arrays[2], arrays[2][end])
    end

    # Formatting prep
    notimesteps = length(storagesoc_comp[1])
    stepspermonth = notimesteps/12

    Xticks = collect(range(stepspermonth/2, notimesteps-(stepspermonth/2), 12))
    Xticklabels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    Xformonthly = collect(range(1, notimesteps, 13))
    Xforweekly = collect(range(1, notimesteps, 52))
    Xfordaily = collect(range(1, notimesteps, 365))

    alpha = [0.4, 0.5]
    linewidthstairs = 3

    # Plotting Data
    for i in [1,2]
    # Row AB
        # Plot A
        # Storage over the entire year with low opacity
        lines!(axa, Xfordaily, storagesoc_comp_daily[i], color=(selectedcolors[i], alpha[i]), linestyle=:dash, label=tslabels[i])
        # Monthly average with high opacity
        stairs!(axa, Xformonthly, storagesoc_comp_monthly[i], step=:post, linewidth=linewidthstairs, color=selectedcolors[i], labels=tslabels[i])

            # Plot B
        # Storage over the entire year with low opacity
        lines!(axb, Xfordaily, storagesoc_geo_daily[i], color=(selectedcolors[i], alpha[i]), linestyle=:dash, label=tslabels[i])
        # Monthly average with high opacity
        stairs!(axb, Xformonthly, storagesoc_geo_monthly[i], step=:post, linewidth=linewidthstairs, color=selectedcolors[i], labels=tslabels[i])

        # Row CD
            # Plot C
        # RES over the entire year with low opacity, weekly
        lines!(axc, Xfordaily, res_comp_daily[i], linestyle=:dash, color=(selectedcolors[i], alpha[i]))
        # RES over the entire year average
        stairs!(axc, Xformonthly, res_comp_monthly[i], step=:post, linewidth=linewidthstairs, color=selectedcolors[i], label=tslabels[i])
        # Revenue over the entire year average
        #stairs!(axcdual, Xformonthly, revenue_comp_monthly[i], step=:post, linestyle=:dash, linewidth=linewidthstairs, color=selectedcolors[i], label=tslabels[i])
            # Plot D
        # RES over the entire year with low opacity, weekly
        lines!(axd, Xfordaily, res_geo_daily[i], linestyle=:dash, color=(selectedcolors[i], alpha[i]))
        # RES over the entire year average
        stairs!(axd, Xformonthly, res_geo_monthly[i], step=:post, linewidth=linewidthstairs, color=selectedcolors[i], label=tslabels[i])
        # Revenue over the entire year average
        #stairs!(axddual, Xformonthly, revenue_geo_monthly[i], step=:post, linestyle=:dot, linewidth=linewidthstairs, color=selectedcolors[i], label=tslabels[i])

    end
    linkxaxes!(axa, axc)
    linkxaxes!(axb, axd)

    linkyaxes!(axc, axd)

    axc.xticks = (Xticks, Xticklabels)
    axd.xticks = (Xticks, Xticklabels)
    axa.xminorticks = Xformonthly
    axb.xminorticks = Xformonthly
    axc.xminorticks = Xformonthly
    axd.xminorticks = Xformonthly
    xlims!(axc, 0, notimesteps)
    xlims!(axd, 0, notimesteps)
    # xlims!(axcdual, 0, notimesteps)
    # xlims!(axddual, 0, notimesteps)

    Legend(fig[3,:], axc, orientation = :horizontal, tellwidth = true, tellheight = true)

    for (label, layout) in zip(["A", "B", "C", "D"], [ga, gb, gc, gd])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    CSV.write("storagelevels.csv", DataFrame([storagesoc_comp_monthly[1], storagesoc_comp_monthly[2], storagesoc_geo_monthly[1], storagesoc_geo_monthly[2], res_comp_monthly[1], res_comp_monthly[2], res_geo_monthly[1], res_geo_monthly[2]], [:HW80, :MW80, :HW80geo, :MW80geo, :HW80RES, :MW80RES, :HW80geoRES, :MW80geoRES]))

    return fig
end

"""
    Durationcurveplot(profiles_filename, ordering) 
    
Imports two csv's with the time series profiles of electricity generators and electricity 
market prices and plots a duration curve of the profiles, generators on the left, prices on 
the right. The DataFrame should have the following columns:

- Price: The price of electricity on the DAM
- Solar: The time series of solar generation
- WindOnshore: The time series of onshore wind generation

The function also needs the ordering to construct the weighting of the profiles. 


"""
function durationcurveplots(datapairs)
    #Importing the data
    fig = Figure()
    slines = nothing
    wlines = nothing
    plines = nothing
    for (i, data) in enumerate(datapairs)
        profiles_filename, ordering = data
        ordering = (CSV.File(ordering, header=1) |> Tables.matrix)
        weights = vec(sum(ordering, dims=1))
        repbyweight(vals,cnts) = [v for (v,c) in zip(vals, cnts) for i in 1:c]


        profiles = CSV.read(profiles_filename, DataFrame)
        profiles[:, :Weights] = round.(Int64,repeat(weights, inner=24))
        profiles = combine(profiles, Not(:Weights) .=> (x -> repbyweight(x, profiles.Weights)) .=> Not(:Weights))

        tail=24*4

        capacityfactors = [sort(profiles[:, gencat], rev=true) for gencat in ["Solar", "WindOnshore"]]
        capacityfactors = [vcat(cf[1:tail], first(denoise(cf[tail:end-tail], factor=10.0)), cf[end-tail:end]) for cf in capacityfactors]
        meancfs = mean.(capacityfactors)
        xmaxcfs = (last.(findmin.([abs.(capacityfactors[i].-meancfs[i]) for i in eachindex(capacityfactors)]))./length(capacityfactors[1]))
        cflabels = ["Solar", "Wind (Onshore)"]
        cfcolor = [:gold, :dodgerblue3]

        prices = sort(profiles[:, :Price], rev=true)
        prices = vcat(prices[1:tail], first(denoise(prices[tail:end-tail], factor=10.0)), prices[end-tail:end])
        meanprice = mean(prices)
        xminprices= last(findmin(abs.(prices.-meanprice)))/length(prices)

        pricelabel = "Price"
        pricecolor = :black

        xticklabels = ["$(round(Int64,i))" for i in range(0, 100, 5)]
        xtickloc = collect(range(1,length(capacityfactors[1]),5))
        yticks = range(0,1,5)
        yticksprice = range(minimum(prices), maximum(prices), 3)
        yticklabelsprice = ["$(round(Int64,i))" for i in range(minimum(prices), maximum(prices), 3)]
        axlim = 1e-2

        ax1 = Axis(
            fig[i,1],
            xlabel="Percentage of total period [%]",
            ylabel="Capacity factor [%]",
            title=profiles_filename[1:2],
            xticks=(xtickloc, xticklabels),
            yticks=(yticks, xticklabels)
        )
        slines = lines!(ax1, capacityfactors[1], labels=cflabels[1], color=cfcolor[1])
        wlines = lines!(ax1, capacityfactors[2], labels=cflabels[2], color=cfcolor[2])
        hlines!(ax1, meancfs, xmax=xmaxcfs, linestyle=:dash, color=cfcolor)

        # Adding legend
        # axislegend(ax1; position=:lb, framevisible=false, labelsize=14)

        ax2 = Axis(
            fig[i,1], 
            xticks=(xtickloc, xticklabels),
            yticks=(yticksprice, yticklabelsprice),
            yaxisposition=:right,
            ylabel="DAM price [€/MWh]"
        )
        hidespines!(ax2)
        hidexdecorations!(ax2)

        # Adding a twin axis for pricedata
        plines = lines!(ax2, prices, label=pricelabel, color=pricecolor)
        hlines!(ax2, mean(prices), xmin=xminprices, linestyle=:dash, color=pricecolor)

        tightlimits!.(ax2)
        tightlimits!.(ax1)
        # axislegend(ax2; position=:rt)
    end
    Legend(fig[3,1], [slines, wlines], ["Solar", "Wind"], halign=:left, orientation = :horizontal, tellwidth = false, tellheight = true)
    Legend(fig[3,1], [plines], ["DAM price"], halign=:right, orientation = :horizontal, tellwidth = false, tellheight = true)
    fig
end


function plotsensitivity(
    results::Vector{DataFrame},
    references
)
    fig = Figure(resolution=(1200,700))

    ga = fig[1,1] = GridLayout()
    gb = fig[1,2] = GridLayout()

    axa = Axis(
        ga[1,1], title="HW80", titlefont=:bold, 
        ylabel="LCOH [€/kg]", xlabel="Deviation from the baseline capital cost", xgridvisible=true, ygridvisible=true
    )
    axb = Axis(
        gb[1,1], title="MW80", titlefont=:bold, xlabel="Deviation from the baseline capital cost", xgridvisible=true, ygridvisible=true
    )
    linkyaxes!(axa, axb)
    hidespines!.([axa, axb])
    hideydecorations!(axb, grid=false)

    # Extracting the data
    components = ["solar", "wind_on", "electrolyzer", "storage", "compressor"]
    complabels = ["Solar", "Wind (onshore)", "Electrolyzer", "Storage", "Compressor"]
    colors = [:gold, :dodgerblue3, :firebrick, :mediumseagreen, :darkorange]
    markers = [:circle, :diamond, :cross, :utriangle, :star4]



    # Iterating over the components
    for (result, reference, axis) in zip(results, references, [axa,axb])
        vlines!(axis, [0], color=:black, linewidth=1, linestyle=:dash)
        hlines!(axis, [reference[:outcome_data]["LCOH"]], color=:black, linewidth=1, linestyle=:dash)

        adjustments = sort!(unique(result[:, :adjustment]))
        xlabels = ["$(round(Int64,100*(adjustment-1)))%" for adjustment in adjustments]

        middleidx = ceil(Int64, length(xlabels)/2)
        xlabels[middleidx] = "Base"

        println(xlabels)
        xticks = [i-(middleidx) for i in range(1,length(xlabels))]
        axis.xticks = (xticks, xlabels)
        for (i, comp) in enumerate(components)
            data = result[result[:, :component] .== comp, :]
            sort!(data, :adjustment)

            LCOH = convert.(Float64, data[:, :LCOH])

            scatterlines!(axis, xticks, LCOH, label=complabels[i], color=colors[i], marker=markers[i], markersize=20, markeralpha=0.5)

        end
    end
    Legend(fig[2,:], axa, orientation = :horizontal, tellwidth = false, tellheight = true)

    return fig
end

## APPENDIX
function makeplot(
    flows::Array{Array{Float64,1},1},
    labels::Array{String,1};
    position::Tuple{Int64,Int64}=(1,1),
    twindata::Union{Array{Float64,1}, Bool}=false,
    twinlabel::String="",
    subsetsize::Union{Int64, Bool}=false
)
    # Removing the flows and labels that are all zeros
    zeros = findall(x->iszero(x), flows)
    parsedflows = deleteat!(copy(flows), zeros)
    parsedlabels = deleteat!(copy(labels), zeros)

    # Taking flow subsets 
    if subsetsize != false
        start = 1#round(Int, length(flows[1])/2)
        parsedflows = map(x->x[start:start+subsetsize-1], parsedflows)
        if isa(twindata, Array)
            twindata = twindata[1:subsetsize]
        end
    end

    fig = Figure()

    ax1 = Axis(fig[position...])
    series!(ax1, parsedflows, labels=parsedlabels)

    # Adding legend
    #axislegend(ax1; position=:lt)

    ax2 = Axis(fig[position...], yaxisposition=:right)
    hidespines!(ax2)
    hidexdecorations!(ax2)

    # Adding a twin axis if needed
    if isa(twindata, Array)
        lines!(ax2, twindata, label=twinlabel, color=:black, linestyle=:dash)
        #axislegend(ax2; position=:rt)
        return fig, ax1, ax2
    end

    Legend(fig[2,1], [ax1.series_list[1]], parsedlabels, halign=:left, orientation = :horizontal, tellwidth = false, tellheight = true)

    return fig, ax1, ax2
end

function powerfig(
    powerflow::Vector{Vector{Float64}},
    powerlabels::Array{String,1};
    twindata::Vector{Float64},
    twinlabel::String="",
    twinaxislabel::String="",
    subsetsize::Union{Int64, Bool}=false
)
    fig, ax1, ax2 = makeplot(powerflow, powerlabels; twindata=twindata, twinlabel=twinlabel, subsetsize=subsetsize)
    ax1.title = "Power flows"
    ax1.ylabel ="Power [MW]"
    ax1.xlabel = "Time [h]"
    ax2.ylabel = twinaxislabel
    return fig
end

function hydrogenfig(
    hydrogenflow::Array{Array{Float64,1},1},
    hydrogenlabels::Array{String,1};
    twindata::Union{Array{Float64,1}, Bool}=false,
    twinlabel::String="",
    twinaxislabel::String="",
    subsetsize::Union{Int64, Bool}=false
)
    fig, ax1, ax2 = makeplot(hydrogenflow, hydrogenlabels; twindata=twindata, twinlabel=twinlabel, subsetsize=subsetsize)
    ax1.title = "Hydrogen flows"
    ax1.ylabel ="Hydrogen massflow [kg]"
    ax1.xlabel = "Time [h]"
    ax2.ylabel = twinaxislabel
    return fig
end

function makeplots(
    data::DataFrame,
    parameters;
    days::Int64=7,
    periods::Int64=52,
    flowsubsets::Dict=Dict()
)
    figs = Dict{Symbol, Any}()

    # Creating a dictionary to map the flow strings to labels
    labeldict = Dict(
        "p_electrolyzer" => "Electrolyzer",
        "p_solar" => "Solar",
        "p_wind_on" => "Wind (onshore)",
        "p_wind_off" => "Wind (offshore)",
        "p_DAM_buy" => "DAM buy",
        "p_DAM_sell" => "DAM sell",
        "p_compressor" => "Compressor",
        "h_electrolyzer" => "Electrolyzer",
        "h_storage_in" => "Storage (in)",
        "h_storage_out" => "Storage (out)",
        "h_storage_soc" => "Storage (SOC)",
        "h_demand_dir" => "DSP feed (direct)",
        "h_demand" => "DSP feed"
    )

    subsetsize = days*24

    if !isempty(flowsubsets)
        powerflownames = flowsubsets[:power]
        hydrogenflownames = flowsubsets[:hydrogen]
    else
        powerflownames = filter(x->startswith(x, "p_"), names(data))
        hydrogenflownames = filter(x->startswith(x, "h_")&&x!="h_storage_soc", names(data))
    end

    # Collecting all the power data
    powerlabels = [labeldict[char] for char in powerflownames]
    powerflow = [data[:, datapoint] for datapoint in powerflownames]
    price_DAM = rep_to_ext(parameters.timeseries.price_DAM, parameters.timeseries.ordering)

    # Collecting all the hydrogen data
    hydrogenlabels = [labeldict[char] for char in hydrogenflownames]
    hydrogenflow = [data[:, datapoint] for datapoint in hydrogenflownames]
    hsoc = data[:, :h_storage_soc]

    # Plotting the flow subset
    figs[:pfig_sub] = powerfig(powerflow, powerlabels, twindata=price_DAM, twinlabel="DAM Price", twinaxislabel="DAM price [€/MWh]", subsetsize=subsetsize)
    figs[:hfig_sub] = hydrogenfig(hydrogenflow, hydrogenlabels; twindata=hsoc, twinlabel="Storage SOC", twinaxislabel="Storage SOC [kg]", subsetsize=subsetsize)

    # Plotting the full flow with a periodic average
    periodic_powerflow = [sum_periodically(pflow, periods=periods) for pflow in powerflow]
    periodic_price_DAM = sum_periodically(price_DAM, periods=periods)
    periodic_hydrogenflow = [sum_periodically(hflow, periods=periods) for hflow in hydrogenflow]
    periodic_hsoc = sum_periodically(data[:, :h_storage_soc])

    figs[:pfig_full] = powerfig(periodic_powerflow, powerlabels, twindata=periodic_price_DAM, twinlabel="DAM Price", twinaxislabel="DAM price [€/MWh]")
    figs[:hfig_full] = hydrogenfig(periodic_hydrogenflow, hydrogenlabels; twindata=periodic_hsoc, twinlabel="Storage SOC", twinaxislabel="Storage SOC [kg]")


    return figs
end
