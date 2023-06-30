using CSV
using DataFrames
using Statistics
using Dates
using Missings

folderpath = joinpath(@__DIR__, "Rawdata/")

# Emission intensity of specific generators [kgCO2/MWh]
EIgenerators = Dict(
    "NL" => Dict(
        "FossilGas" => 393.24,
        "FossilHardCoal" => 1183.07
    ),
    "ES" => Dict(
        "FossilGas" => 390.43,
        "FossilHardCoal" => 1044.87,
        "FossilOil" => 840
    ) 
)

Region = "NL"
TimeZone = "CET"

DAfiles = [
    "DAMp_20_21.csv", 
    "DAMp_21_22.csv",
    "DAMp_22_23.csv"
]
AGfiles = [
    "AGpPT_20_21.csv", 
    "AGpPT_21_22.csv",
    "AGpPT_22_23.csv"
    ]

columnnamedict = Dict(
    "MTU" => :Time,
    "MTU (CET/CEST)" => :Time,
    "Day-ahead Price [EUR/MWh]" => :DAprice,
    "Solar  - Actual Aggregated [MW]" => :Solar,
    "Wind Onshore  - Actual Aggregated [MW]" => :WindOnshore,
    "Wind Offshore  - Actual Aggregated [MW]" => :WindOffshore,
    "Biomass  - Actual Aggregated [MW]" => :Biomass,
    "Geothermal  - Actual Aggregated [MW]" => :Geothermal,
    "Fossil Gas  - Actual Aggregated [MW]" => :FossilGas,
    "Fossil Hard coal  - Actual Aggregated [MW]" => :FossilHardCoal,
    "Fossil Oil  - Actual Aggregated [MW]" => :FossilOil,
    "Hydro Pumped Storage  - Actual Aggregated [MW]" => :HydroPumpedStorageGeneration,
    "Hydro Pumped Storage  - Actual Consumption [MW]" => :HydroPumpedStorageConsumption,
    "Hydro Run-of-river and poundage  - Actual Aggregated [MW]" => :HydroRunOfRiver,
    "Hydro Water Reservoir  - Actual Aggregated [MW]" => :HydroWaterReservoir,
    "Marine  - Actual Aggregated [MW]" => :Marine,
    "Nuclear  - Actual Aggregated [MW]" => :Nuclear,
    "Other  - Actual Aggregated [MW]" => :Other,
    "Other renewable  - Actual Aggregated [MW]" => :OtherRenewable,
    "Waste  - Actual Aggregated [MW]" => :Waste
)

function fillgapswithmean(
    df::DataFrame,
    columnnames::Array{String,1}
)
    for columnname in columnnames
        #println("Working on column $(columnname)")
        array = df[!, columnname]
        for i in 1:length(df[:, columnname])
            if ismissing(array[i])
                #println("found missing at row $i: $(array[i])")
                if i == 1
                    #println("Replacing missing value at row $(i) with mean")
                    array[i] = mean(skipmissing(array))
                #if both neighbors are there, we can take the mean
                elseif sum(ismissing.(array[[i-1, i+1]])) == 0 
                    #println("Replacing missing value at row $(i) with mean of surrounding values")
                    array[i] = mean(array[[i-1, i+1]])
                elseif i <= 24
                    #println("Replacing missing value at row $(i) with mean")
                    array[i] = mean(skipmissing(array))
                elseif !(ismissing(df[i-24, columnname]))
                    array[i] = df[i-24, columnname]
                else
                    throw("No neighbors or previous day")
                end
            end
        end
    end
    return df
end

function timecolumntodatetime!(
    df::DataFrame,
    timecolumn::String
)
    datetimestrings = [DateTime.(split(row, " - ")[1], "dd.mm.yyyy HH:MM") for row in df[:, timecolumn]]
    df[!, timecolumn] = datetimestrings
end

function preparedataframe(
    df::DataFrame,
    timecolumn::String,
    datacolumns::AbstractVector,
    normalize::Bool
)
    if isempty(datacolumns)
        datacolumns = filter((name->name!=timecolumn), names(df))
        datacolumns = [colname for colname in names(df) if (any(val->isa(val, Float64), df[:,colname]) && any(val->val>0, skipmissing(df[:,colname])))]
    end
    df = df[:, [timecolumn, datacolumns...]]
    display(df)
    fillgapswithmean(df, datacolumns)
    disallowmissing!(df)
    timecolumntodatetime!(df, timecolumn)
    firstquarterly = isquarterhour(df, timecolumn)
    if !isnothing(firstquarterly)
        if !isnothing(firstquarterly)
            println("Splitting quarterly")
            hourlydf = df[1:firstquarterly-1, :]
            quarterlydf = df[firstquarterly:end, :]
            quarterlydf = meanrowintervals(quarterlydf, 4, timecolumn, datacolumns)
            df = vcat(hourlydf, quarterlydf)
        end
    end
    if normalize
        df = normalizecolumns(df, [timecolumn])
    end
    rename!(df, getindex.(Ref(columnnamedict), names(df)))
    return df
end

function dataframefromcsvs(
    csvfilenames::Array{String,1},
    timecolumn::String;
    normalize::Bool=false,
    datacolumns::AbstractVector=[]
)
    dfs = [
        CSV.File(
            joinpath(folderpath, csvfilename),
            missingstring=["N/A", "", " ", "n/e"],
            typemap=Dict(Int=>Float64)
        ) |> DataFrame
        for csvfilename in csvfilenames
    ]
    #println("Total timesteps: ", [if isquarterhour(df) size(df,1)/4 else size(df,1) end for df in dfs])
    dfs = [preparedataframe(df, timecolumn, datacolumns, normalize) for df in dfs]
    return vcat(dfs...)
end

function findrowsofmissingvalues(
    df::DataFrame
)
    for columnname in names(df)
        #println("Finding missing values in column $(columnname)")
        for i in 1:length(df[:, columnname])
            if ismissing(df[i, columnname])
                #println("Missing value at row $(i)")
            elseif isa(df[i, columnname], String)
                #println("String at row $(i)")
            elseif isnan(df[i, columnname])
                #println("NaN value at row $(i)")
            end
        end
    end
end

function normalizecolumns(
    df::DataFrame,
    except::Union{Array{String, 1},Array{Symbol, 1}} = []
)   
    for columnname in names(select(df, Not(except)))
        #println("Normalizing column $(columnname), with mean value: $(mean(df[:, columnname]))")
        maxval = maximum(df[:, columnname])
        #println("Max value: $(maxval)")
        if maxval != 0
            df[:, columnname] = df[:, columnname] ./ maxval
        end
    end
    return df
end

function meanrowintervals(
    df::DataFrame, 
    interval::Int,
    timecolumn::String,
    datacolumns::Vector{String}
)
    df.group = repeat(1:(nrow(df)/interval), inner = interval)
    grouped_df = groupby(df, :group)
    mean_df = combine(grouped_df, timecolumn .=> first, datacolumns .=> mean, renamecols=false)
    return select(mean_df, Not(:group))
end

function isquarterhour(df, timecol)
    tatrail = df[:, timecol][1:end-1]
    talead = df[:, timecol][2:end]
    deltat = (talead .- tatrail)
    quarterly = deltat .== Dates.Minute(15)

    # # Dealing with summer/wintertime transitions
    # forwardskip = findall((deltat .== Dates.Minute(120)) .|| (deltat .==  Dates.Minute(75)))
    # backwardskip = findall((deltat .==  -Dates.Minute(45)) .|| (deltat .==  Dates.Minute(0)))
    # if length(forwardskip) > 1 || length(backwardskip) > 1
    #     throw("Unexplained time differences: $(df[:, timecol][forwardskip]) and $(df[:, timecol][backwardskip])")
    # elseif (length(forwardskip) == length(backwardskip) == 1) || (length(forwardskip) == length(backwardskip) == 0)
    #     println("Found transition at $(df[:, timecol][forwardskip]) and $(df[:, timecol][backwardskip]), no changes")
    # elseif length(forwardskip) == 1 && length(backwardskip) == 0
    #     println("Found transition at $(df[:, timecol][forwardskip]), deleting")
    #     if quarterly[forwardskip]
    #         datetimes = df[collect((forwardskip-3):forwardskip), timecol]
    #         throw(datetimes)
    #         insert!(df.columns, collect(collect((forwardskip-3):forwardskip), timecol))
    #     else
    #         deleteat!(df , idx)
    #     end
    # elseif length(forwardskip) == 0 && length(backwardskip) == 1
    #     println("Found transition at $(df[:, timecol][backwardskip]), deleting")
    #     if quarterly[backwardskip]
    #         deleteat!(df , collect(idx:(idx+3)))
    #     else
    #         deleteat!(df , idx)
    #     end
    # end
    return findfirst(quarterly)
end

function savedataframetocsv(
    df::DataFrame, 
    csvfilename::String
)
    CSV.write(csvfilename, df)
end

function getEFdf(
    df::DataFrame,
    region::String
)
    consumptioncolumns = filter(name->occursin("Consumption", name), names(df))
    datacolumns = filter((name->name!="Time"), names(df))
    totalgendata = sum(Matrix(df[:,datacolumns]),dims=2)
    println(describe(df))
    if !isempty(consumptioncolumns)
        totalconsdata = sum(Matrix(df[:,consumptioncolumns]),dims=2)
    else
        totalconsdata = 0
    end
    df[:, "TotalGeneration"] .= totalgendata .- totalconsdata
    df[:, "TotalEmissions"] = sum([df[:, genkey] .* EIgenerators[region][genkey] for genkey in keys(EIgenerators[region])])
    df[:, "EI"] = df[:, "TotalEmissions"] ./ df[:, "TotalGeneration"]
    println(describe(df))
    select!(df, ["Time", "EI"])
    return df
end

function createEFcsv(
    generationfilenames::Array{String,1},
    region::String;
    timecol::String="MTU"
)
    gendf = dataframefromcsvs(
        generationfilenames, 
        timecol
    )
    EFdf = getEFdf(gendf, region)
    return EFdf
end

function createtimeseriesdf(
    dayaheadfilenames::Array{String,1},
    generationfilenames::Array{String,1},
    region::String;
    outfilename::String=""
)
    println("Fetching day-ahead data")
    dadf = dataframefromcsvs(
        dayaheadfilenames, 
        "MTU (CET/CEST)",
        datacolumns=["Day-ahead Price [EUR/MWh]"]
    )
    println("Fetching capacity factor data")
    gendf = dataframefromcsvs(
        generationfilenames, 
        "MTU",
        datacolumns=["Solar  - Actual Aggregated [MW]", 
            "Wind Onshore  - Actual Aggregated [MW]", 
            "Wind Offshore  - Actual Aggregated [MW]"
        ],
        normalize=true
    )
    println("Fetching emissions data")
    EIdf = getEFdf(
        dataframefromcsvs(
            generationfilenames, 
            "MTU"),
        region
    )
    for (tda, tgen, tei) in zip(dadf[:, :Time], gendf[:, :Time], EIdf[:, :Time])
        if !(tda == tgen == tei)
            throw("Times don't match at: $(tda) vs $(tgen) vs $(tei)")
        end
    end
    fulldf = hcat(dadf, select(gendf, Not(:Time)), select(EIdf, Not(:Time)))
    return fulldf
end

function createcsvfile()
    df = createtimeseriesdf(
        Region.*"_".*DAfiles, 
        Region.*"_".*AGfiles, 
        Region)
    CSV.write("$(Region)timeseries.csv", df)
end

# Function for checking a problem in the hourly daily rep_to_ext translation equivalence
function findprob()
    sumvalue = 0

    sprob = "0.01 os_electrolyzer[1]*c_electrolyzer + 0.01 os_electrolyzer[2]*c_electrolyzer + 0.01 os_electrolyzer[3]*c_electrolyzer + 0.01 os_electrolyzer[4]*c_electrolyzer + 0.01 os_electrolyzer[5]*c_electrolyzer + 0.01 os_electrolyzer[6]*c_electrolyzer + 0.01 os_electrolyzer[7]*c_electrolyzer + 0.01 os_electrolyzer[8]*c_electrolyzer + 0.01 os_electrolyzer[9]*c_electrolyzer + 0.01 os_electrolyzer[10]*c_electrolyzer + 0.01 os_electrolyzer[11]*c_electrolyzer + 0.01 os_electrolyzer[12]*c_electrolyzer + 0.01 os_electrolyzer[13]*c_electrolyzer + 0.01 os_electrolyzer[14]*c_electrolyzer + 0.01 os_electrolyzer[15]*c_electrolyzer + 0.01 os_electrolyzer[16]*c_electrolyzer + 0.01 os_electrolyzer[17]*c_electrolyzer + 0.01 os_electrolyzer[18]*c_electrolyzer + 0.01 os_electrolyzer[19]*c_electrolyzer + 0.01 os_electrolyzer[20]*c_electrolyzer + 0.01 os_electrolyzer[21]*c_electrolyzer + 0.01 os_electrolyzer[22]*c_electrolyzer + 0.01 os_electrolyzer[23]*c_electrolyzer + 0.01 os_electrolyzer[24]*c_electrolyzer + p_DAM_buy[1] - p_DAM_sell[1] - p_compressor[1] - 0.24000000000000007 c_electrolyzer + p_DAM_buy[2] - p_DAM_sell[2] - p_compressor[2] + p_DAM_buy[3] - p_DAM_sell[3] - p_compressor[3] + p_DAM_buy[4] - p_DAM_sell[4] - p_compressor[4] + p_DAM_buy[5] - p_DAM_sell[5] - p_compressor[5] + p_DAM_buy[6] - p_DAM_sell[6] - p_compressor[6] + p_DAM_buy[7] - p_DAM_sell[7] - p_compressor[7] + p_DAM_buy[8] - p_DAM_sell[8] - p_compressor[8] + p_DAM_buy[9] - p_DAM_sell[9] - p_compressor[9] + p_DAM_buy[10] - p_DAM_sell[10] - p_compressor[10] + p_DAM_buy[11] - p_DAM_sell[11] - p_compressor[11] + p_DAM_buy[12] - p_DAM_sell[12] - p_compressor[12] + p_DAM_buy[13] - p_DAM_sell[13] - p_compressor[13] + p_DAM_buy[14] - p_DAM_sell[14] - p_compressor[14] + p_DAM_buy[15] - p_DAM_sell[15] - p_compressor[15] + p_DAM_buy[16] - p_DAM_sell[16] - p_compressor[16] + p_DAM_buy[17] - p_DAM_sell[17] - p_compressor[17] + p_DAM_buy[18] - p_DAM_sell[18] - p_compressor[18] + p_DAM_buy[19] - p_DAM_sell[19] - p_compressor[19] + p_DAM_buy[20] - p_DAM_sell[20] - p_compressor[20] + p_DAM_buy[21] - p_DAM_sell[21] - p_compressor[21] + p_DAM_buy[22] - p_DAM_sell[22] - p_compressor[22] + p_DAM_buy[23] - p_DAM_sell[23] - p_compressor[23] + p_DAM_buy[24] - p_DAM_sell[24] - p_compressor[24], 0.01 os_electrolyzer[1]*c_electrolyzer + 0.01 os_electrolyzer[2]*c_electrolyzer + 0.01 os_electrolyzer[3]*c_electrolyzer + 0.01 os_electrolyzer[4]*c_electrolyzer + 0.01 os_electrolyzer[5]*c_electrolyzer + 0.01 os_electrolyzer[6]*c_electrolyzer + 0.01 os_electrolyzer[7]*c_electrolyzer + 0.01 os_electrolyzer[8]*c_electrolyzer + 0.01 os_electrolyzer[9]*c_electrolyzer + 0.01 os_electrolyzer[10]*c_electrolyzer + 0.01 os_electrolyzer[11]*c_electrolyzer + 0.01 os_electrolyzer[12]*c_electrolyzer + 0.01 os_electrolyzer[13]*c_electrolyzer + 0.01 os_electrolyzer[14]*c_electrolyzer + 0.01 os_electrolyzer[15]*c_electrolyzer + 0.01 os_electrolyzer[16]*c_electrolyzer + 0.01 os_electrolyzer[17]*c_electrolyzer + 0.01 os_electrolyzer[18]*c_electrolyzer + 0.01 os_electrolyzer[19]*c_electrolyzer + 0.01 os_electrolyzer[20]*c_electrolyzer + 0.01 os_electrolyzer[21]*c_electrolyzer + 0.01 os_electrolyzer[22]*c_electrolyzer + 0.01 os_electrolyzer[23]*c_electrolyzer + 0.01 os_electrolyzer[24]*c_electrolyzer + p_DAM_buy[1] - p_DAM_sell[1] - p_compressor[1] - 0.24000000000000007 c_electrolyzer + p_DAM_buy[2] - p_DAM_sell[2] - p_compressor[2] + p_DAM_buy[3] - p_DAM_sell[3] - p_compressor[3] + p_DAM_buy[4] - p_DAM_sell[4] - p_compressor[4] + p_DAM_buy[5] - p_DAM_sell[5] - p_compressor[5] + p_DAM_buy[6] - p_DAM_sell[6] - p_compressor[6] + p_DAM_buy[7] - p_DAM_sell[7] - p_compressor[7] + p_DAM_buy[8] - p_DAM_sell[8] - p_compressor[8] + p_DAM_buy[9] - p_DAM_sell[9] - p_compressor[9] + p_DAM_buy[10] - p_DAM_sell[10] - p_compressor[10] + p_DAM_buy[11] - p_DAM_sell[11] - p_compressor[11] + p_DAM_buy[12] - p_DAM_sell[12] - p_compressor[12] + p_DAM_buy[13] - p_DAM_sell[13] - p_compressor[13] + p_DAM_buy[14] - p_DAM_sell[14] - p_compressor[14] + p_DAM_buy[15] - p_DAM_sell[15] - p_compressor[15] + p_DAM_buy[16] - p_DAM_sell[16] - p_compressor[16] + p_DAM_buy[17] - p_DAM_sell[17] - p_compressor[17] + p_DAM_buy[18] - p_DAM_sell[18] - p_compressor[18] + p_DAM_buy[19] - p_DAM_sell[19] - p_compressor[19] + p_DAM_buy[20] - p_DAM_sell[20] - p_compressor[20] + p_DAM_buy[21] - p_DAM_sell[21] - p_compressor[21] + p_DAM_buy[22] - p_DAM_sell[22] - p_compressor[22] + p_DAM_buy[23] - p_DAM_sell[23] - p_compressor[23] + p_DAM_buy[24] - p_DAM_sell[24] - p_compressor[24]"
    swork = "0.01 os_electrolyzer[1]*c_electrolyzer + p_DAM_buy[1] - p_DAM_sell[1] - p_compressor[1] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[2]*c_electrolyzer + p_DAM_buy[2] - p_DAM_sell[2] - p_compressor[2] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[3]*c_electrolyzer + p_DAM_buy[3] - p_DAM_sell[3] - p_compressor[3] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[4]*c_electrolyzer + p_DAM_buy[4] - p_DAM_sell[4] - p_compressor[4] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[5]*c_electrolyzer + p_DAM_buy[5] - p_DAM_sell[5] - p_compressor[5] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[6]*c_electrolyzer + p_DAM_buy[6] - p_DAM_sell[6] - p_compressor[6] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[7]*c_electrolyzer + p_DAM_buy[7] - p_DAM_sell[7] - p_compressor[7] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[8]*c_electrolyzer + p_DAM_buy[8] - p_DAM_sell[8] - p_compressor[8] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[9]*c_electrolyzer + p_DAM_buy[9] - p_DAM_sell[9] - p_compressor[9] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[10]*c_electrolyzer + p_DAM_buy[10] - p_DAM_sell[10] - p_compressor[10] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[11]*c_electrolyzer + p_DAM_buy[11] - p_DAM_sell[11] - p_compressor[11] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[12]*c_electrolyzer + p_DAM_buy[12] - p_DAM_sell[12] - p_compressor[12] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[13]*c_electrolyzer + p_DAM_buy[13] - p_DAM_sell[13] - p_compressor[13] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[14]*c_electrolyzer + p_DAM_buy[14] - p_DAM_sell[14] - p_compressor[14] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[15]*c_electrolyzer + p_DAM_buy[15] - p_DAM_sell[15] - p_compressor[15] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[16]*c_electrolyzer + p_DAM_buy[16] - p_DAM_sell[16] - p_compressor[16] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[17]*c_electrolyzer + p_DAM_buy[17] - p_DAM_sell[17] - p_compressor[17] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[18]*c_electrolyzer + p_DAM_buy[18] - p_DAM_sell[18] - p_compressor[18] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[19]*c_electrolyzer + p_DAM_buy[19] - p_DAM_sell[19] - p_compressor[19] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[20]*c_electrolyzer + p_DAM_buy[20] - p_DAM_sell[20] - p_compressor[20] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[21]*c_electrolyzer + p_DAM_buy[21] - p_DAM_sell[21] - p_compressor[21] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[22]*c_electrolyzer + p_DAM_buy[22] - p_DAM_sell[22] - p_compressor[22] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[23]*c_electrolyzer + p_DAM_buy[23] - p_DAM_sell[23] - p_compressor[23] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[24]*c_electrolyzer + p_DAM_buy[24] - p_DAM_sell[24] - p_compressor[24] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[1]*c_electrolyzer + p_DAM_buy[1] - p_DAM_sell[1] - p_compressor[1] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[2]*c_electrolyzer + p_DAM_buy[2] - p_DAM_sell[2] - p_compressor[2] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[3]*c_electrolyzer + p_DAM_buy[3] - p_DAM_sell[3] - p_compressor[3] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[4]*c_electrolyzer + p_DAM_buy[4] - p_DAM_sell[4] - p_compressor[4] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[5]*c_electrolyzer + p_DAM_buy[5] - p_DAM_sell[5] - p_compressor[5] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[6]*c_electrolyzer + p_DAM_buy[6] - p_DAM_sell[6] - p_compressor[6] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[7]*c_electrolyzer + p_DAM_buy[7] - p_DAM_sell[7] - p_compressor[7] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[8]*c_electrolyzer + p_DAM_buy[8] - p_DAM_sell[8] - p_compressor[8] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[9]*c_electrolyzer + p_DAM_buy[9] - p_DAM_sell[9] - p_compressor[9] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[10]*c_electrolyzer + p_DAM_buy[10] - p_DAM_sell[10] - p_compressor[10] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[11]*c_electrolyzer + p_DAM_buy[11] - p_DAM_sell[11] - p_compressor[11] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[12]*c_electrolyzer + p_DAM_buy[12] - p_DAM_sell[12] - p_compressor[12] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[13]*c_electrolyzer + p_DAM_buy[13] - p_DAM_sell[13] - p_compressor[13] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[14]*c_electrolyzer + p_DAM_buy[14] - p_DAM_sell[14] - p_compressor[14] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[15]*c_electrolyzer + p_DAM_buy[15] - p_DAM_sell[15] - p_compressor[15] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[16]*c_electrolyzer + p_DAM_buy[16] - p_DAM_sell[16] - p_compressor[16] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[17]*c_electrolyzer + p_DAM_buy[17] - p_DAM_sell[17] - p_compressor[17] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[18]*c_electrolyzer + p_DAM_buy[18] - p_DAM_sell[18] - p_compressor[18] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[19]*c_electrolyzer + p_DAM_buy[19] - p_DAM_sell[19] - p_compressor[19] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[20]*c_electrolyzer + p_DAM_buy[20] - p_DAM_sell[20] - p_compressor[20] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[21]*c_electrolyzer + p_DAM_buy[21] - p_DAM_sell[21] - p_compressor[21] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[22]*c_electrolyzer + p_DAM_buy[22] - p_DAM_sell[22] - p_compressor[22] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[23]*c_electrolyzer + p_DAM_buy[23] - p_DAM_sell[23] - p_compressor[23] - 0.01 c_electrolyzer, 0.01 os_electrolyzer[24]*c_electrolyzer + p_DAM_buy[24] - p_DAM_sell[24] - p_compressor[24] - 0.01 c_electrolyzer"
    
    probarray = strip.(split(sprob, (',' , '+', '-')))
    workarray = strip.(split(swork, (',' , '+', '-')))
    for substr in workarray
        try
            deleteat!(probarray, findfirst(x -> x == substr, probarray))
        catch
            if substr=="0.01 c_electrolyzer"
                sumvalue += 0.01
            end
            println("$substr not found")
        end
    end
    println(probarray)
    println(sumvalue)
end