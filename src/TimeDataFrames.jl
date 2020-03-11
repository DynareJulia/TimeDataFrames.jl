module TimeDataFrames

using CSV
using DataFrames
using Periods

export TimeDataFrame

struct TimeDataFrame
    data::DataFrame
    periods::Vector{Period}
    continuous::Bool
    frequency::Frequency
end;

TimeDataFrame() = TimeDataFrame(DataFrame(), Vector{Period}(), true, Undated)

TimeDataFrame(frequency::Frequency) = TimeDataFrame(DataFrame(), Vector{Period}(), true, frequency)

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

names(df::TimeDataFrame) = names(getfield(df, :data))

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
        minperiod = (period1[1] < period2[1]) ? periods1[1] : periods2[1]
        maxperiod = (period1[n1] > period2[n2]) ? periods1[n1] : periods2[n2]
        # Missing expanding number of periods in data1 if necessary
    end
#    if !(col_ind in DataFrames.names(data1))
#    end
#    setindex!(data1, data2[!, 1], periods2[1]:periods2[n2], col_ind)
    TimeDataFrame(data1, periods2, continuous2, frequency2) 
end

Base.getindex(df::TimeDataFrame, id1::Integer, id2::Integer) = getindex(getfield(df, :data), id1, id2)
Base.getindex(df::TimeDataFrame, idx::CartesianIndex{2}) = df[idx[1], idx[2]]
Base.view(df::TimeDataFrame, idx::CartesianIndex{2}) = view(df, idx[1], idx[2])
Base.setindex!(df::TimeDataFrame, val, idx::CartesianIndex{2}) =
    (df[idx[1], idx[2]] = val)

Base.broadcastable(tdf::TimeDataFrame) = tdf
struct TimeDataFrameStyle <: Base.Broadcast.BroadcastStyle end

Base.Broadcast.BroadcastStyle(::Type{<:TimeDataFrame}) =
    TimeDataFrameStyle()

Base.Broadcast.BroadcastStyle(::TimeDataFrameStyle, ::Base.Broadcast.BroadcastStyle) = TimeDataFrameStyle()
Base.Broadcast.BroadcastStyle(::Base.Broadcast.BroadcastStyle, ::TimeDataFrameStyle) = TimeDataFrameStyle()
Base.Broadcast.BroadcastStyle(::TimeDataFrameStyle, ::TimeDataFrameStyle) = TimeDataFrameStyle()

function Base.copy(bc::Base.Broadcast.Broadcasted{TimeDataFrameStyle})
    ndim = length(axes(bc))
    if ndim != 2
        throw(DimensionMismatch("cannot broadcast a time data frame into $ndim dimensions"))
    end

    data 
    bcf = Base.Broadcast.flatten(bc)
    colnames = unique!([_names(df) for df in bcf.args if df isa AbstractDataFrame])
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
    df = DataFrame()
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
        df[!, colnames[1][i]] = col
    end
    return df
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
