module TimeDataFrames

using CSV
using DataFrames
using ExtendedDates

import Base: eachcol, eachrow, show

export TimeDataFrame, periods, firstperiod, lastperiod, dataframe, innerjoin, outerjoin, lag, lead, ncol, nrow

mutable struct TimeDataFrame
    data::DataFrame
    periods::AbstractVector{ExtendedDates.SimpleDate}
    continuous::Bool
    function TimeDataFrame(data::DataFrame,
                           periods::AbstractVector{T},
                           continuous::Bool) where T <: ExtendedDates.SimpleDate
        new(data, periods, continuous)
    end
end;

    
TimeDataFrame() = TimeDataFrame(DataFrame(), Vector{T}(), true) where T <: ExtendedDates.SimpleDate

function TimeDataFrame(dataframe::AbstractDataFrame, firstperiod::T;
                       copycols::Bool=true) where T <: ExtendedDates.SimpleDate
    periods = range(firstperiod, length=DataFrames.nrow(dataframe), step = typeof(firstperiod - firstperiod)(1))
    TimeDataFrame(DataFrame(dataframe; copycols), periods, true) 
end

function TimeDataFrame(filename::String, firstperiod::T) where T <: ExtendedDates.SimpleDate
    data = DataFrame(CSV.File(filename))
    continuous = true
    periods = range(firstperiod, length=size(data, 1), step = typeof(firstperiod - firstperiod)(1))
    TimeDataFrame(data, periods, true) 
end

TimeDataFrame(firstperiod::T, pairs::Pair{Symbol,<:Any}...;
              makeunique::Bool=false, copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(pairs; makeunique, copycols),
                                firstperiod)

TimeDataFrame(firstperiod::T, pairs::Pair{<:AbstractString,<:Any}...;
              makeunique::Bool=false, copycols::Bool=true) where T <: ExtendedDates.SimpleDate =  
                  TimeDataFrame(DataFrame(pairs; makeunique, copycols),
                                firstperiod)

# these two are needed as a workaround Tables.jl dispatch
TimeDataFrame(pairs::AbstractVector{<:Pair}, firstperiod::T; makeunique::Bool=false,
              copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(pairs..., makeunique=makeunique, copycols=copycols),
                                firstperiod)

TimeDataFrame(pairs::NTuple{N, Pair}, firstperiod::T; makeunique::Bool=false,
              copycols::Bool=true) where {N,T <: ExtendedDates.SimpleDate} =
                  TimeDataFrame(DataFrame(pairs..., makeunique=makeunique, copycols=copycols),
                                firstperiod) 

TimeDataFrame(d::AbstractDict, firstperiod::T; copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
    TimeDataFrame(DataFrame(d; copycols), firstperiod)

#=
TimeDataFrame(frequency::ExtendedDates.Frequency, firstperiod::ExtendedDates.SimpleDate; kwargs...)
    if isempty(kwargs)
        TimeDataFrame()
    else
        cnames = Symbol[]
        columns = Any[]
        copycols = true
        for (kw, val) in kwargs
            if kw == :copycols
                if val isa Bool
                    copycols = val
                else
                    throw(ArgumentError("the `copycols` keyword argument must be Boolean"))
                end
            else
                push!(cnames, kw)
                push!(columns, val)
            end
        end
        TimeDataFrame(DataFrame(columns, Index(cnames), copycols=copycols), Undated, 1)
    end
end
=#

TimeDataFrame(columns::AbstractVector, cnames::AbstractVector{Symbol},
              firstperiod::T; makeunique::Bool=false, copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(columns, cnames; makeunique, copycols), firstperiod)

TimeDataFrame(columns::AbstractVector, cnames::AbstractVector{<:AbstractString},
              firstperiod::T; makeunique::Bool=false, copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(columns, Symbol.(cnames), makeunique=makeunique, copycols=copycols),
                                firstperiod)

TimeDataFrame(columns::AbstractVector{<:AbstractVector}, firstperiod::T,
              cnames::AbstractVector{Symbol}=gennames(length(columns));
              makeunique::Bool=false, copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(collect(AbstractVector, columns),
                                          Index(convert(Vector{Symbol}, cnames), makeunique=makeunique),
                                          copycols=copycols), firstperiod)

TimeDataFrame(columns::AbstractVector{<:AbstractVector},
              cnames::AbstractVector{<:AbstractString}, firstperiod::T;
              makeunique::Bool=false, copycols::Bool=true) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(columns, Symbol.(cnames); makeunique=makeunique, copycols=copycols),
                                firstperiod)

TimeDataFrame(columns::NTuple{N, AbstractVector}, cnames::NTuple{N, Symbol}, firstperiod::T;
              makeunique::Bool=false, copycols::Bool=true) where {N, T <: ExtendedDates.SimpleDate} =
                  TimeDataFrame(DataFrame(collect(AbstractVector, columns), collect(Symbol, cnames),
                                          makeunique=makeunique, copycols=copycols), firstperiod)

TimeDataFrame(columns::NTuple{N, AbstractVector}, cnames::NTuple{N, AbstractString}, firstperiod::T;
          makeunique::Bool=false, copycols::Bool=true) where {N, T <: ExtendedDates.SimpleDate} =
              TimeDataFrame(DataFrame(columns, Symbol.(cnames); makeunique=makeunique, copycols=copycols),
                            firstperiod)

TimeDataFrame(columns::NTuple{N, AbstractVector}, firstperiod::T; copycols::Bool=true) where {N, T <: ExtendedDates.SimpleDate} =
    TimeDataFrame(DataFrame(collect(AbstractVector, columns), gennames(length(columns)),
                            copycols=copycols), firstperiod)

TimeDataFrame(columns::AbstractMatrix, firstperiod::T,
          cnames::AbstractVector{Symbol} = gennames(size(columns, 2)); makeunique::Bool=false) where T <: ExtendedDates.SimpleDate =
              TimeDataFrame(DataFrame(AbstractVector[columns[:, i] for i in 1:size(columns, 2)], cnames,
                                      makeunique=makeunique, copycols=false), firstperiod)

TimeDataFrame(columns::AbstractMatrix, cnames::AbstractVector{<:AbstractString}, firstperiod::T;
              makeunique::Bool=false) where T <: ExtendedDates.SimpleDate =
                  TimeDataFrame(DataFrame(columns, Symbol.(cnames); makeunique=makeunique), firstperiod)

TimeDataFrame(column_eltypes::AbstractVector{T}, cnames::AbstractVector{Symbol},
              firstperiod::P, nrows::Integer=0; makeunique::Bool=false) where {T<:Type, P <: ExtendedDates.SimpleDate} =
                  TimeDataFrame(DataFrame(column_eltypes, cnames, nrows; makeunique), firstperiod)

TimeDataFrame(column_eltypes::AbstractVector{T}, cnames::AbstractVector{S}, firstperiod::P,
              nrows::Integer=0; makeunique::Bool=false) where {T <: Type, P <: ExtendedDates.SimpleDate, S <: AbstractString} =
                  TimeDataFrame(DataFrame(column_eltypes, Symbol.(cnames), nrows; makeunique=makeunique),
                                firstperiod)

# Redefining getproperty breaks tdf.data !
function Base.getproperty(tdf::TimeDataFrame, symbol::Symbol)
    data = getfield(tdf, :data)
    x = getproperty(data, symbol)
    periods = getfield(tdf, :periods)
    continuous = getfield(tdf, :continuous)
    df = DataFrame()
    df[!, symbol] = x
    TimeDataFrame(df,  periods, continuous)
end

Base.names(tdf::TimeDataFrame) = DataFrames.names(DataFrames.index(getfield(tdf, :data)))
_names(tdf::TimeDataFrame) = DataFrames._names(DataFrames.index(getfield(tdf, :data)))

function Base.setproperty!(tdf::TimeDataFrame, symbol::Symbol, x::AbstractVector)
    data = getfield(tdf, :data)
    periods = getfield(tdf, :periods)
    data[!, symbol] = x
    continuous = getfield(tdf, :continuous)
    TimeDataFrame(data, periods, continuous)
end

function Base.setproperty!(tdf1::TimeDataFrame, col_ind::Symbol, tdf2::TimeDataFrame)
    n1 = nrow(tdf1)
    n2 = nrow(tdf2)
    continuous2 = getfield(tdf2, :continuous)
    data1 = getfield(tdf1, :data)
    data2 = getfield(tdf2, :data)
    periods2 = getfield(tdf2, :periods)
    if ncol(tdf1) == 0
        insertcols!(data1, 1, col_ind => missings(n2))
    else
        continuous1 = getfield(tdf1, :continuous)
        periods1 = getfield(tdf1, :periods)
        if typeof(periods1) != typeof(periods2)
            error("The frequency must be the same in both TimeDataFrame")
        end
        minperiod = (periods1[1] < periods2[1]) ? periods1[1] : periods2[1]
        maxperiod = (periods1[n1] > periods2[n2]) ? periods1[n1] : periods2[n2]
        # Missing expanding number of periods in data1 if necessary
    end
#    if !(col_ind in DataFrames.names(data1))
#    end
    setindex!(data1, data2[!, 1], ! , col_ind)
#    rename!(data1, [col_ind])
    setfield!(tdf1, :data, data1)
    setfield!(tdf1, :periods, periods2)
    setfield!(tdf1, :continuous, continuous2)
end

Base.getindex(df::TimeDataFrame, id1::Integer, id2::Integer) = getindex(getfield(df, :data), id1, id2)
Base.getindex(df::TimeDataFrame, idx::CartesianIndex{2}) = df[idx[1], idx[2]]
Base.getindex(df::TimeDataFrame, ::typeof(!), col_ind::Symbol) = getindex(getfield(df, :data), !, col_ind)
Base.getindex(df::TimeDataFrame, ::typeof(!), col_ind::Union{Signed, Unsigned}) = getindex(getfield(df, :data), !, col_ind)

Base.view(df::TimeDataFrame, idx::CartesianIndex{2}) = view(df, idx[1], idx[2])
Base.setindex!(tdf::TimeDataFrame, val, idx::CartesianIndex{2}) =
    (tdf[idx[1], idx[2]] = val)
Base.setindex!(tdf::TimeDataFrame, x::AbstractArray{Float64, 1}, ::typeof(!), col_ind::Symbol) =
               setindex!(getfield(tdf, :data), x, !, col_ind)
Base.setindex!(tdf::TimeDataFrame, x::AbstractArray{Union{Missing, Float64}, 1}, ::typeof(!), col_ind::Symbol) =
               setindex!(getfield(tdf, :data), x, !, col_ind)
Base.setindex!(tdf::TimeDataFrame, x::AbstractArray{Float64, 1}, ::typeof(!), col_ind::Union{Signed, Unsigned}) =
    setindex!(getfield(tdf, :data), x, !, col_ind)

Base.broadcastable(tdf::TimeDataFrame) = tdf
struct TimeDataFrameStyle <: Base.Broadcast.BroadcastStyle end

Base.Broadcast.BroadcastStyle(::Type{<:TimeDataFrame}) =
    TimeDataFrameStyle()

Base.Broadcast.BroadcastStyle(::TimeDataFrameStyle, ::Base.Broadcast.BroadcastStyle) = TimeDataFrameStyle()
Base.Broadcast.BroadcastStyle(::Base.Broadcast.BroadcastStyle, ::TimeDataFrameStyle) = TimeDataFrameStyle()
Base.Broadcast.BroadcastStyle(::TimeDataFrameStyle, ::TimeDataFrameStyle) = TimeDataFrameStyle()

function copyto_widen!(res::AbstractVector{T}, bc::Base.Broadcast.Broadcasted,
                       pos, col) where T
    for i in pos:length(axes(bc)[1])
        val = bc[CartesianIndex(i, col)]
        S = typeof(val)
        if S <: T || promote_type(S, T) <: T
            res[i] = val
        else
            newres = similar(Vector{promote_type(S, T)}, length(res))
            copyto!(newres, 1, res, 1, i-1)
            newres[i] = val
            return copyto_widen!(newres, bc, i + 1, col)
        end
    end
    return res
end

function getcolbc(bcf::Base.Broadcast.Broadcasted{Style}, colind) where {Style}
    # we assume that bcf is already flattened and unaliased
    newargs = map(bcf.args) do x
        Base.Broadcast.extrude(x isa TimeDataFrame ? x[!, colind] : x)
    end
    Base.Broadcast.Broadcasted{Style}(bcf.f, newargs, bcf.axes)
end

function Base.copy(bc::Base.Broadcast.Broadcasted{TimeDataFrameStyle})
    ndim = length(axes(bc))
    if ndim != 2
        throw(DimensionMismatch("cannot broadcast a time data frame into $ndim dimensions"))
    end

    
    bcf = Base.Broadcast.flatten(bc)
    first_tdf = true
    local periods, continuous
    colnames = []
    for tdf in bcf.args
        if tdf isa TimeDataFrame
            if first_tdf
                periods = getfield(tdf, :periods)
                continuous = getfield(tdf, :continuous)
                first_tdf = false
            elseif getfield(tdf, :periods) != periods
                error("Time data frames don't have the same periods")
            elseif getfield(tdf, :continuous) != continuous
                error("TimeDataFrames don't have the same continuous status")
            end
        end
    end            
    colnames = unique!([_names(df) for df in bcf.args if df isa TimeDataFrame])
    if length(colnames) != 1
        wrongnames = setdiff(union(colnames...), intersect(colnames...))
        if isempty(wrongnames)
            throw(ArgumentError("Column names in broadcasted data frames " *
                                "must have the same order"))
        else
            msg = join(wrongnames, ", ", " and ")
            throw(ArgumentError("Column names in broadcasted data frames must match. " *
                                "Non matching column names are $msg"))
        end
    end
    nrows = length(axes(bcf)[1])
    tdf = TimeDataFrame(DataFrame(), periods, continuous)
    for i in axes(bcf)[2]
        if nrows == 0
            col = Any[]
        else
            bcf′ = getcolbc(bcf, i)
            v1 = bcf′[CartesianIndex(1, i)]
            startcol = similar(Vector{typeof(v1)}, nrows)
            startcol[1] = v1
            col = copyto_widen!(startcol, bcf′, 2, i)
        end
        tdf[!, colnames[1][i]] = col
    end
    return tdf
end




Base.ndims(::TimeDataFrame) = 2
Base.ndims(::Type{<:TimeDataFrame}) = 2
index(df::TimeDataFrame) = getfield(getfield(df, :data), :colindex)
_columns(df::TimeDataFrame) = getfield(getfield(df, :data), :columns)

# note: these type assertions are required to pass tests
ncol(df::TimeDataFrame) = length(index(df))
nrow(df::TimeDataFrame) = ncol(df) > 0 ? length(_columns(df)[1])::Int : 0
import Base.size
Base.size(df::TimeDataFrame) = (nrow(df), ncol(df))
function Base.size(df::TimeDataFrame, i::Integer)
    if i == 1
        nrow(df)
    elseif i == 2
        ncol(df)
    else
        throw(ArgumentError("DataFrames only have two dimensions"))
    end
end

function innerjoin(d1::TimeDataFrame, d2::TimeDataFrame)
    continuous1 = getfield(d1, :continuous)
    continuous2 = getfield(d2, :continuous)
    data1 = getfield(d1, :data)
    data2 = getfield(d2, :data)
    periods1 = getfield(d1, :periods)
    periods2 = getfield(d2, :periods)
    if typeof(periods1) != typeof(periods2)
        error("innerjoin: both TimeDataFrames must have the same frequency")
    end
    return TimeDataFrame(sort!(DataFrames.innerjoin(data1, data2, on=:Column1),1), intersect(periods1, periods2), true)
end

function outerjoin(d1::TimeDataFrame, d2::TimeDataFrame)
    continuous1 = getfield(d1, :continuous)
    continuous2 = getfield(d2, :continuous)
    data1 = getfield(d1, :data)
    data2 = getfield(d2, :data)
    periods1 = getfield(d1, :periods)
    periods2 = getfield(d2, :periods)
    if typeof(periods1) != typeof(periods2)
        error("innerjoin: both TimeDataFrames must have the same frequency")
    end
    return TimeDataFrame(sort!(DataFrames.outerjoin(data1, data2, on=:Column1),1), union(periods1, periods2), true)
end

function Base.show(io::IO, tdf::TimeDataFrame)
    df = getfield(tdf, :data)
    dfcopy = copy(df)
    periods = getfield(tdf, :periods)
    insertcols!(dfcopy, 1, :Periods => periods)
    show(io, dfcopy, show_row_number = false, eltypes = false, summary = false)
end

function Base.show(tdf::TimeDataFrame)
    df = getfield(tdf, :data)
    dfcopy = copy(df)
    periods = getfield(tdf, :periods)
    insertcols!(dfcopy, 1, :Periods => periods)
    show(dfcopy, show_row_number = false, eltypes = false, summary = false)
end

function Base.isequal(tdf1::TimeDataFrame, tdf2::TimeDataFrame)
    isequal(getfield(tdf1, :data), getfield(tdf2, :data)) || return false
    isequal(getfield(tdf1, :periods), getfield(tdf2, :periods)) || return false
    getfield(tdf1, :continuous) == getfield(tdf2, :continuous) || return false
    return true
end

include("accessors.jl")
include("dataframe_functions.jl")
include("timeseries_functions.jl")

end # module
