module TimeDataFrames

using CSV
using DataFrames
using Periods

struct TimeDataFrame
    data::DataFrame
    continuous::Bool
    frequency::Periods.Frequency
end;

function TimeDataFrame(dataframe::AbstractDataFrame, frequency::Periods.Frequency, firstperiod)
    data = dataframe
    continuous = true
    data.Periods = [Periods(firstperiod + i - 1, 0, Periods.Year) for i in 1:dataframe.nrow] 
    TimeDataFrame(dataframe, true, frequency) 
end

function TimeDataFrame(filename::String, frequency::Periods.Frequency, firstperiod)
    data = DataFrame(CSV.File(filename))
    continuous = true
    data.periods = [Periods.SimplePeriod(firstperiod + i - 1, 0, Periods.Year) for i in 1:nrow(data)] 
    TimeDataFrame(data, true, frequency) 
end

end # module
