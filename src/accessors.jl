"""
    periods(tdf::TimeDataFrame)
returns the periods of the TimeDataFrame
"""
periods(tdf::TimeDataFrame) = getfield(tdf, :periods)

"""
    firstperiod(tdf::TimeDataFrame)
returns the first period of the TimeDataFrame
"""
firstperiod(tdf::TimeDataFrame) = periods(tdf)[1]

"""
    lasttperiod(tdf::TimeDataFrame)
returns the last period of the TimeDataFrame
"""
lastperiod(tdf::TimeDataFrame) = periods(tdf)[end]

"""
    dataframe(tdf::TimeDataFrame)
returns the DataFrame inside the TimeDataFrame
"""
dataframe(tdf::TimeDataFrame) = getfield(tdf, :data)

"""
    iscontinuous(tdf::TimeDataFrame)
returns whether the TimeDataFrame is continuous
"""
iscontinuous(tdf::TimeDataFrame) = getfield(tdf, :continuous)
