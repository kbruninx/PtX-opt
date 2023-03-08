using CSV
using DataFrames
using Statistics

filename = "Actual Generation per Production Type_201901010000-202001010000.csv"

divisordict = Dict("Solar  - Actual Aggregated [MW]" => 4608, "Wind Onshore  - Actual Aggregated [MW]" => 3436, "Wind Offshore  - Actual Aggregated [MW]" => 957)

function fillgapswithmean(
    df::DataFrame, 
    columnname::String
)
    for i in 1:length(df[:, columnname])
        if isnothing(df[i, columnname])
            println("Replacing missing value at row $(i) with mean of surrounding values")
            df[i, columnname] = mean(df[[i-1, i+1], columnname])
        end
    end
    return df
end

function dataframefromcsv(
    csvfilename::String
)
    df = CSV.read(csvfilename, DataFrame, types = Float64)
    return df
end

function normalizecolumns(
    df::DataFrame
)   
    for columnname in names(df)
        maxval = maximum(df[:, columnname])
        df[:, columnname] = df[:, columnname] ./ maxval
    end
    return df
end

function meanrowintervals(
    df::DataFrame, 
    interval::Int
)
    df.group = repeat(1:8761, inner = interval)
    grouped_df = groupby(df, :group)
    mean_df = combine(grouped_df, valuecols(grouped_df) .=> mean)
    return select(mean_df, Not(:group))
end

function savedataframetocsv(
    df::DataFrame, 
    csvfilename::String
)
    CSV.write(csvfilename, df)
end

function main_csvhandler(
    infilename::String, 
    outfilename::String,
    interval::Int
)
    df = dataframefromcsv(infilename)
    df = normalizecolumns(df)
    df = meanrowintervals(df, interval)
    savedataframetocsv(df, outfilename)
end