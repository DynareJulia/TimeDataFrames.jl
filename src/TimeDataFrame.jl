module TimeDataFrame

import DataFrame
import Periods

struct TimeDataFrame
    data::DataFrame
    continuous::Bool
    frequency::Periods.Frequency
end;

function TimeDataFrame(dataframe::DataFrame, frequency::Periods.Frequency, firstperiod)
    data = dataframe
    continuous = True
    data.Periods = [Periods(firstperido + i - 1, 0, Periods.Year) for i in 1:dataframe.nrow] 
    TimeDataFrame(dataframe, true, frequency) 
end

end # module
