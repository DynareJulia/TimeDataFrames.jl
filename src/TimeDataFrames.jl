module TimeDataFrames

using CSV
using DataFrames
using Periods

export TimeDataFrame

struct TimeDataFrame
    data::DataFrame
    continuous::Bool
    frequency::Frequency
end;

TimeDataFrame() = TimeDataFrame(DataFrame(), true, Undated)

TimeDataFrame(frequency::Frequency) = TimeDataFrame(DataFrame(), true, frequency)

function TimeDataFrame(dataframe::AbstractDataFrame, frequency::Frequency, firstperiod)
    data = dataframe
    continuous = true
    data.Periods = [Period(firstperiod + i - 1, 0, Year) for i in 1:dataframe.nrow] 
    TimeDataFrame(dataframe, true, frequency) 
end

function TimeDataFrame(filename::String, frequency::Frequency, firstperiod)
    data = DataFrame(CSV.File(filename))
    continuous = true
    data.periods = [Period(firstperiod + i - 1, 0, Year) for i in 1:nrow(data)] 
    TimeDataFrame(data, true, frequency) 
end

# Redefining getproperty breaks tdf.data !
function Base.getproperty(tdf::TimeDataFrame, symbol::Symbol)
    data = getfield(tdf, :data)
    x = getproperty(data, symbol)
    periods = getproperty(data, :periods)
    continuous = getfield(tdf, :continuous)
    frequency = getfield(tdf, :frequency)
    TimeDataFrame(DataFrame(periods = periods, symbol = x ),  continuous, frequency)
end

Base.setproperty!(tdf::TimeDataFrame, symbol::Symbol, x::AbstractArray) = setproperty!(getfield(tdf, :data), symbol, x)

end # module
