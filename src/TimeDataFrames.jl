module TimeDataFrames

using CSV
using DataFrames
using Periods

export TimeDataFrame

mutable struct TimeDataFrame
    data::DataFrame
    periods::Vector{Period}
    continuous::Bool
    frequency::Frequency
end;

TimeDataFrame() = TimeDataFrame(DataFrame(), Vector{Period}(), true,
                                Undated)

TimeDataFrame(frequency::Frequency) = TimeDataFrame(DataFrame(),
                                                    Vector{Period}(), true, frequency)

function TimeDataFrame(dataframe::AbstractDataFrame, frequency::Frequency, firstperiod)
    periods = [Period(firstperiod + i - 1, 0, frequency)  for i in 1:DataFrames.nrow(dataframe)]
    TimeDataFrame(dataframe, periods, true, frequency) 
end

function TimeDataFrame(filename::String, frequency::Frequency, firstperiod)
    data = DataFrame(CSV.File(filename))
    continuous = true
    periods = [Period(firstperiod + i - 1, 0, Year) for i in 1:size(data, 1)] 
    TimeDataFrame(data, periods, true, frequency) 
end

TimeDataFrame(df::DataFrame, frequency::Frequency, firstperiod;
              copycols::Bool=true) =
                  TimeDataFrame(DataFrame(df; copycols),
                                frequency, firstperiod)

TimeDataFrame(frequency::Frequency, firstperiod, pairs::Pair{Symbol,<:Any}...;
              makeunique::Bool=false, copycols::Bool=true) =
                  TimeDataFrame(DataFrame(pairs; makeunique, copycols),
                                frequency, firstperiod)

TimeDataFrame(frequency::Frequency, firstperiod, pairs::Pair{<:AbstractString,<:Any}...;
              makeunique::Bool=false, copycols::Bool=true) =
                  TimeDataFrame(DataFrame(pairs; makeunique, copycols),
                                frequency, firstperiod)

# these two are needed as a workaround Tables.jl dispatch
TimeDataFrame(pairs::AbstractVector{<:Pair}, frequency::Frequency, firstperiod; makeunique::Bool=false,
              copycols::Bool=true) =
                  TimeDataFrame(DataFrame(pairs..., makeunique=makeunique, copycols=copycols),
                                frequency, firstperiod)

TimeDataFrame(pairs::NTuple{N, Pair}, frequency::Frequency, firstperiod; makeunique::Bool=false,
              copycols::Bool=true) where {N} =
                  TimeDataFrame(DataFrame(pairs..., makeunique=makeunique, copycols=copycols),
                                frequency, firstperiod) 

TimeDataFrame(d::AbstractDict, frequency::Frequency, firstperiod; copycols::Bool=true) =
    TimeDataFrame(DataFrame(d; copycols), frequency, firstperiod)

#=
TimeDataFrame(frequency::Frequency, firstperiod; kwargs...)
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

TimeDataFrame(columns::AbstractVector, cnames::AbstractVector{Symbol}, frequency::Frequency,
              firstperiod; makeunique::Bool=false, copycols::Bool=true) =
                  TimeDataFrame(DataFrame(columns, cnames; makeunique, copycols), frequency, firstperiod)

TimeDataFrame(columns::AbstractVector, cnames::AbstractVector{<:AbstractString}, frequency::Frequency,
              firstperiod; makeunique::Bool=false, copycols::Bool=true) =
                  TimeDataFrame(DataFrame(columns, Symbol.(cnames), makeunique=makeunique, copycols=copycols),
                                frequency, firstperiod)

TimeDataFrame(columns::AbstractVector{<:AbstractVector}, frequency::Frequency, firstperiod,
              cnames::AbstractVector{Symbol}=gennames(length(columns));
              makeunique::Bool=false, copycols::Bool=true) =
                  TimeDataFrame(DataFrame(collect(AbstractVector, columns),
                                          Index(convert(Vector{Symbol}, cnames), makeunique=makeunique),
                                          copycols=copycols), frequency, firstperiod)

TimeDataFrame(columns::AbstractVector{<:AbstractVector},
              cnames::AbstractVector{<:AbstractString}, frequency::Frequency, firstperiod;
              makeunique::Bool=false, copycols::Bool=true) =
                  TimeDataFrame(DataFrame(columns, Symbol.(cnames); makeunique=makeunique, copycols=copycols),
                                frequency, firstperiod)

TimeDataFrame(columns::NTuple{N, AbstractVector}, cnames::NTuple{N, Symbol}, frequency::Frequency, firstperiod;
              makeunique::Bool=false, copycols::Bool=true) where {N} =
                  TimeDataFrame(DataFrame(collect(AbstractVector, columns), collect(Symbol, cnames),
                                          makeunique=makeunique, copycols=copycols), frequency, firstperiod)

TimeDataFrame(columns::NTuple{N, AbstractVector}, cnames::NTuple{N, AbstractString}, frequency::Frequency, firstperiod;
          makeunique::Bool=false, copycols::Bool=true) where {N} =
              TimeDataFrame(DataFrame(columns, Symbol.(cnames); makeunique=makeunique, copycols=copycols), frequency,
                            firstperiod)

TimeDataFrame(columns::NTuple{N, AbstractVector}, frequency::Frequency, firstperiod; copycols::Bool=true) where {N} =
    TimeDataFrame(DataFrame(collect(AbstractVector, columns), gennames(length(columns)),
                            copycols=copycols), frequency, firstperiod)

TimeDataFrame(columns::AbstractMatrix, frequency::Frequency, firstperiod,
          cnames::AbstractVector{Symbol} = gennames(size(columns, 2)); makeunique::Bool=false) =
              TimeDataFrame(DataFrame(AbstractVector[columns[:, i] for i in 1:size(columns, 2)], cnames,
                                      makeunique=makeunique, copycols=false), frequency, firstperiod)

TimeDataFrame(columns::AbstractMatrix, cnames::AbstractVector{<:AbstractString};
          makeunique::Bool=false) =
              TimeDataFrame(DataFrame(columns, Symbol.(cnames); makeunique=makeunique), frequency, firstperiod)

TimeDataFrame(column_eltypes::AbstractVector{T}, cnames::AbstractVector{Symbol},
              frequency::Frequency, firstperiod, nrows::Integer=0; makeunique::Bool=false) where T<:Type =
                  TimeDataFrame(DataFrame(column_eltypes, cnames, nrows; makeunique), frequency, firstperiod)

TimeDataFrame(column_eltypes::AbstractVector{<:Type},
              cnames::AbstractVector{<:AbstractString},
              frequency::Frequency, firstperiod, nrows::Integer=0; makeunique::Bool=false) =
                  TimeDataFrame(DataFrame(column_eltypes, Symbol.(cnames), nrows; makeunique=makeunique), frequency,
                                firstperiod)

# Redefining getproperty breaks tdf.data !
function Base.getproperty(tdf::TimeDataFrame, symbol::Symbol)
    data = getfield(tdf, :data)
    x = getproperty(data, symbol)
    periods = getfield(tdf, :periods)
    continuous = getfield(tdf, :continuous)
    frequency = getfield(tdf, :frequency)
    df = DataFrame()
    df[!, symbol] = x
    TimeDataFrame(df,  periods, continuous, frequency)
end

Base.names(tdf::TimeDataFrame) = DataFrames.names(DataFrames.index(getfield(tdf, :data)))
_names(tdf::TimeDataFrame) = DataFrames._names(DataFrames.index(getfield(tdf, :data)))

function Base.setproperty!(tdf::TimeDataFrame, symbol::Symbol, x::AbstractVector)
    data = getfield(tdf, :data)
    periods = getfield(tdf, :periods)
    data[!, symbol] = x
    continuous = getfield(tdf, :continuous)
    frequency = getfield(tdf, :frequency)
    TimeDataFrame(data, periods, continuous, frequency)
end

function Base.setproperty!(tdf1::TimeDataFrame, col_ind::Symbol, tdf2::TimeDataFrame)
    n1 = nrow(tdf1)
    n2 = nrow(tdf2)
    frequency2 = getfield(tdf2, :frequency)
    continuous2 = getfield(tdf2, :continuous)
    data1 = getfield(tdf1, :data)
    data2 = getfield(tdf2, :data)
    periods2 = getfield(tdf2, :periods)
    if ncol(tdf1) == 0
        insertcols!(data1, 1, col_ind => missings(n2))
    else
        frequency1 = getfield(tdf1, :frequency)
        continuous1 = getfield(tdf1, :continuous)
        if frequency1 != frequency2
            error("The frequency must be the same in both TimeDataFrame")
        end
        periods1 = getfield(tdf1, :periods)
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
    setfield!(tdf1, :frequency, frequency2)
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
    local periods, continuous, frequency
    colnames = []
    for tdf in bcf.args
        if tdf isa TimeDataFrame
            if first_tdf
                periods = getfield(tdf, :periods)
                continuous = getfield(tdf, :continuous)
                frequency = getfield(tdf, :frequency)
                first_tdf = false
            elseif getfield(tdf, :periods) != periods
                error("Time data frames don't have the same periods")
            elseif getfield(tdf, :continuous) != continuous
                error("TimeDataFrames don't have the same continuous status")
            elseif getfield(tdf, :frequency) != frequency
                error("TimdeDataFrames don't have the same frequency")
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
    tdf = TimeDataFrame(DataFrame(), periods, continuous, frequency)
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

end # module
