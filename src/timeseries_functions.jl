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

function addmissingfirst!(tdf, k::DatePeriod)
    kv = k.value
    tmp = Vector{Union{Float64, Missing}}(undef, size(tdf,1) + kv)
    for c in eachcol(tdf)
        tmp .= append!(c, repeat([missing], kv))
        c .= circshift(tmp, kv)
    end
end

function addmissingend!(tdf, k::DatePeriod)
    kv = k.value
    for c in eachcol(tdf)
        append!(c, repeat([missing], kv))
    end
end

function addperiodsfirst!(periods, k::DatePeriod)
    kv = k.value
    append!(periods,
            collect(range(periods[1] - typeof(k)(1), length = kv, step=-typeof(k)(1))))
    tmp = circshift(periods, kv)
    periods .= tmp
end

addperiodsend!(periods, k::DatePeriod) =
    append!(periods,
            collect(range(periods[end] + typeof(k)(1), length = k.value, step=typeof(k)(1))))

function align!(tdf, tdfs...)
    periods = TimeDataFrames.periods(tdf)
    p1 = periods[1]
    p2 = periods[end]
    tp = typeof(periods[1])
    for t in tdfs
        periods1 = TimeDataFrames.periods(t)
        if typeof(periods1[1]) != tp
            error("TimeDataFrames must have the same frequency")
        end
        p1 = min(p1, periods1[1])
        p2 = max(p2, periods1[end])
    end
    zeroperiod = periods[1] - periods[1]
    for t in push!([tdfs...], tdf)
        periods1 = TimeDataFrames.periods(t)
        m = periods1[1] - p1
        if m > zeroperiod
            allowmissing!(dataframe(t))
            addmissingfirst!(t, m)
            addperiodsfirst!(periods1, m)
        end
        m = p2 - periods1[end]
        if m > zeroperiod
            allowmissing!(dataframe(t))
            addmissingend!(t, m)
            addperiodsend!(periods1, m)
        end
    end
end
