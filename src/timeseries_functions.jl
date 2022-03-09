function TimeDataFrame(f, tdf::TimeDataFrame, args...)
    df = dataframe(tdf)
    f1(x) = f(x, args)
    return TimeDataFrame(map(f, eachcol(df)), names(df), firstperiod(tdf))
end

function lag(x::AbstractVector{T}, k::Int64) where T
    n = length(x)
    y = Vector{Union{T, Missing}}(undef, n)
    view(y, 1:k) .= missing
    view(y, k+1:n) .= view(x, 1:n-k)
    return y
end


lag(tdf::TimeDataFrame, k::Int64) = TimeDataFrame(lag, tdf, k)
lag(x) = lag(x, 1)

function lead(x::AbstractVector{T}, k::Int64) where T
    n = length(x)
    y = Vector{Union{T, Missing}}(undef, n)
    view(y, 1:n-k) .= view(x, k+1:n)
    view(y, n-k+1:n) .= missing
    return y
end

lead(tdf::TimeDataFrame, k::Int64) = TimeDataFrame(lead, tdf, k)
lead(x) = lead(x, 1)

