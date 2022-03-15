Base.eachcol(tdf::TimeDataFrame) = eachcol(dataframe(tdf))
Base.eachrow(tdf::TimeDataFrame) = eachrow(dataframe(tdf))
ncol(tdf::TimeDataFrame) = DataFrames.ncol(dataframe(tdf))
nrow(tdf::TimeDataFrame) = DataFrames.nrow(dataframe(tdf))
