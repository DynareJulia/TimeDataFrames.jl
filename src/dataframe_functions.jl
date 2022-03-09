Base.eachcol(tdf::TimeDataFrames) = eachcol(dataframe(tdf))
Base.eachrow(tdf::TimeDataFrames) = eachrow(dataframe(tdf))
ncol(tdf::TimeDataFrames) = ncol(dataframe(tdf))
nrow(tdf::TimeDataFrames) = nrow(dataframe(tdf))
