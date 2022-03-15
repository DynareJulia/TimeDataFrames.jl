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
<<<<<<< HEAD

"""
    continuous(tdf::TimeDataFrame)
returns whether the TimeDataFrame is continuous
"""
continuous(tdf::TimeDataFrame) = getfield(tdf, :continuous)

=======
"""
    continuous(tdf::TimeDataFrame)
returns whether the TimeDataFrame is time continuous
"""
continuous(tdf::TimeDataFrame) = getfield(tdf, :continuous)
>>>>>>> 6c1e2a7 (update)
