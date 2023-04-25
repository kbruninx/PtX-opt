using CSV
using DataFrames
using Statistics
using Dates
using Missings

DAfiles = [
    "Day-ahead Prices_202001010000-202101010000.csv", 
    "Day-ahead Prices_202101010000-202201010000.csv",
    "Day-ahead Prices_202201010000-202301010000.csv"
]
AGfiles = [
    "Actual Generation per Production Type_202001010000-202101010000.csv", 
    "Actual Generation per Production Type_202101010000-202201010000.csv",
    "Actual Generation per Production Type_202201010000-202301010000.csv"
    ]

columnnamedict = Dict(
    "MTU" => :Time,
    "MTU (CET/CEST)" => :Time,
    "Day-ahead Price [EUR/MWh]" => :DAprice,
    "Solar  - Actual Aggregated [MW]" => :Solar,
    "Wind Onshore  - Actual Aggregated [MW]" => :WindOnshore,
    "Wind Offshore  - Actual Aggregated [MW]" => :WindOffshore
)

function fillgapswithmean(
    df::DataFrame,
    columnnames::Array{String,1}
)
    println(describe(df))
    for columnname in columnnames
        array = df[!, columnname]
        for i in 1:length(df[:, columnname])
            if ismissing(array[i])
                println("found missing at row $i")
                neighbors = array[[i-1, i+1]]
                #if both neighbors are there, we can take the mean
                if sum(ismissing.(neighbors)) == 0 
                    println("Replacing missing value at row $(i) with mean of surrounding values")
                    array[i] = mean(neighbors)
                elseif !(ismissing(df[i-24, columnname]))
                    array[i] = df[i-24, columnname]
                else
                    throw("No neighbors or previous day")
                end
            end
        end
    end
    println(describe(df))
    return df
end

function timecolumntodatetime(
    df::DataFrame,
    timecolumn::String
)
    datetimestrings = [DateTime(split(row, " - ")[1], "dd.mm.yyyy HH:MM") for row in df[:, timecolumn]]
    df[!, timecolumn] = datetimestrings
end

function preparedataframe(
    df::DataFrame,
    timecolumn::String,
    datacolumns::Array{String,1},
    normalize::Bool
)
    println(first(df, 8))
    fillgapswithmean(df, datacolumns)
    disallowmissing!(df)
    if isquarterhour(df)
        println("Quarterly resolution for:")
        println(first(df, 8))
        df = meanrowintervals(df, 4, timecolumn, datacolumns)
        println("Changed to: ")
        println(first(df, 8))
    end
    if normalize
        df = normalizecolumns(df, [timecolumn])
    end
    timecolumntodatetime(df, timecolumn)
    rename!(df, getindex.(Ref(columnnamedict), names(df)))
    return df
end

function dataframefromcsvs(
    csvfilenames::Array{String,1},
    timecolumn::String,
    datacolumns::Array{String,1};
    normalize::Bool=false
)
    dfs = [
        CSV.File(
            csvfilename,
            select=[timecolumn, datacolumns...],
            missingstring=["N/A", "", " "],
            typemap=Dict(Int=>Float64)
        ) |> DataFrame
        for csvfilename in csvfilenames
    ]
    dfs = [preparedataframe(df, timecolumn, datacolumns, normalize) for df in dfs]
    return vcat(dfs...)
end

function findrowsofmissingvalues(
    df::DataFrame
)
    for columnname in names(df)
        println("Finding missing values in column $(columnname)")
        for i in 1:length(df[:, columnname])
            if ismissing(df[i, columnname])
                println("Missing value at row $(i)")
            elseif isa(df[i, columnname], String)
                println("String at row $(i)")
            elseif isnan(df[i, columnname])
                println("NaN value at row $(i)")
            end
        end
    end
end

function normalizecolumns(
    df::DataFrame,
    except::Union{Array{String, 1},Array{Symbol, 1}} = []
)   
    for columnname in names(select(df, Not(except)))
        println("Normalizing column $(columnname), with mean value: $(mean(df[:, columnname]))")
        maxval = maximum(df[:, columnname])
        println("Max value: $(maxval)")
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

function isquarterhour(df)
    if nrow(df)>9000 && (nrow(df)%4)==0
        return true
    else 
        return false
    end
end

function savedataframetocsv(
    df::DataFrame, 
    csvfilename::String
)
    CSV.write(csvfilename, df)
end

function createtimeseriescsv(
    dayaheadfilenames::Array{String,1},
    generationfilenames::Array{String,1},
    outfilename::String
)
    dadf = dataframefromcsvs(
        dayaheadfilenames, 
        "MTU (CET/CEST)",
        ["Day-ahead Price [EUR/MWh]"]
    )
    gendf = dataframefromcsvs(
        generationfilenames, 
        "MTU",
        ["Solar  - Actual Aggregated [MW]", 
            "Wind Onshore  - Actual Aggregated [MW]", 
            "Wind Offshore  - Actual Aggregated [MW]"
        ],
        normalize=true
    )
    fulldf = innerjoin(
        dadf, gendf,
        on = :Time
    )
    return fulldf
end