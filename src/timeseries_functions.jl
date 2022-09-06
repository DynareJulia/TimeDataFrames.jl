function TimeDataFrame(f::Function, tdf::TimeDataFrame, args...)
    df = dataframe(tdf)
    f1(x) = f(x, args...)
    fp = firstperiod(tdf)
    return TimeDataFrame(mapcols(f1, df), fp)
end

function lag(x::AbstractVector{T}, k::Int64) where T
    n = length(x)
    y = Vector{Union{T, Missing}}(undef, n)
    view(y, 1:k) .= missing
    view(y, k+1:n) .= view(x, 1:n-k)
    return y
end


lag(tdf::TimeDataFrame, k::Int64) = (k > 0) ? TimeDataFrame(lag, tdf, k) : TimeDataFrame(lead, tdf, -k)
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

function iscontinuous(tdfs...)
    for tdf in tdfs
        iscontinuous(tdf) || return false
    end
    return true
end

function iscontinuous(tdf, tdfs...)
    iscontinuous(tdf) || return false
    for t in tdfs
        iscontinuous(t) || return false
    end
    return true
end

function havesameperiods(tdfs::TimeDataFrame...)
    havesameperiods(Val(iscontinuous(tdfs)), tdfs)
end

function havesameperiods(::Val{true}, tdf::TimeDataFrame, tdfs::TimeDataFrame...)
    p1 = firstperiod(tdf)
    p2 = lastperiod(tdf)
    for t in tdfs
        firstperiod(t) == p1 || return false
        lastperiod(t) == p2 || return false
    end
    return true
end

function havesameperiods(::Val{false}, tdf::TimeDataFrame, tdfs::TimeDataFrame...)
    p = periods(tdf)
    for t in tdfs
        periods(t) == p || return false
    end
    return true
end

function align_!(tdf, tdfs...)
    (iscontinuous(tdf, tdfs...)) ? align_!(Val(true), tdf, tdfs...) : align_!(Val(false), tdf, tdfs...)
end

function align_!(::Val{true}, tdf, tdfs...)
    p = periods(tdf)
    p1 = p[1]
    p2 = p[end]
    tp = typeof(p[1])
    for t in tdfs
        p = periods(t)
        if typeof(p[1]) != tp
            error("TimeDataFrames must have the same frequency")
        end
        p1 = min(p1, p[1])
        p2 = max(p2, p[end])
    end
    zeroperiod = p[1] - p[1]
    for t in push!([tdfs...], tdf)
        p = periods(t)
        m = p[1] - p1
        if m > zeroperiod
            allowmissing!(dataframe(t))
            addmissingfirst!(t, m)
            addperiodsfirst!(p, m)
        end
        m = p2 - p[end]
        if m > zeroperiod
            allowmissing!(dataframe(t))
            addmissingend!(t, m)
            addperiodsend!(p, m)
        end
    end
end

function align_!(::Val{false}, tdf, tdfs...)
    error("Not implemented")
end

function align(tdf, tdfs...)
    tdf1 = copy(tdf)
    tdfs1 = [copy(t) for t in tdfs]
    align_!(tdf1, tdfs1...)
    return(tdf1, tdfs1)
end

function align(tdf1::TimeDataFrame, tdf2::TimeDataFrame)
    tdf1a = copy(tdf1)
    tdf2a = copy(tdf2)
    align_!(tdf1a, tdf2a)
    return (tdf1a, tdf2a)
end

function align!(tdf, tdfs...)
    align_!(tdf, tdfs...)
    return(tdf, tdfs)
end

import Base: +, -
    
function +(x1::TimeDataFrame, y1::TimeDataFrame)
    (TimeDataFrames.ncol(x1) != TimeDataFrames.ncol(y1)) && error("TimeDataFrames must have the same number of columns")
    (x, y) = align(x1, y1)
    if TimeDataFrames.ncol(x) == 1
        return TimeDataFrame( Matrix(x) + Matrix(y), ["$(names(x)[1])_+_$(names(y)[1])"], firstperiod(x))
    else
        names(x) != names(y) && error("Multicolumn TimeDataFrames must have the same column names")
        S = Matrix{Union{Float64, Missing}}(undef, size(x))
        dfx = dataframe(x)
        dfy = dataframe(y)
        for (i, n) in enumerate(names(x))
            view(S, :, i) = getproperty(dfx, Symbol(n)) .+ getproperty(dfy, Symbol(n))
        end
        return TimeDataFrame(S, names(x), firstperiod(x))
    end
end

function -(x1::TimeDataFrame, y1::TimeDataFrame)
    (TimeDataFrames.ncol(x1) != TimeDataFrames.ncol(y1)) && error("TimeDataFrames must have the same number of columns")
    (x, y) = align(x1, y1)
    if TimeDataFrames.ncol(x) == 1
        return TimeDataFrame( Matrix(x) - Matrix(y), ["$(names(x)[1])_-_$(names(y)[1])"], firstperiod(x))
    else
        names(x) != names(y) && error("Multicolumn TimeDataFrames must have the same column names")
        S = Matrix{Union{Float64, Missing}}(undef, size(x))
        dfx = dataframe(x)
        dfy = dataframe(y)
        for (i, n) in enumerate(names(x))
            view(S, :, i) = getproperty(dfx, Symbol(n)) .- getproperty(dfy, Symbol(n))
        end
        return TimeDataFrame(S, names(x), firstperiod(x))
    end
end

function index(tdf::TimeDataFrame, base::ExtendedDates.SimpleDate)
    baseindex = findall(base .== periods(tdf))
    return TimeDataFrame(x -> 100 .* x./x[baseindex], tdf)
end
                        
