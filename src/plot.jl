using Plots
import Plots: plot

function plot(tdf::TimeDataFrame; label = [], title = "")
    x = [string(p) for p in periods(tdf)]
    if TimeDataFrames.ncol(tdf) > 1 && isempty(label)
        label = names(tdf)
    end
    for (i,c) in enumerate(eachcol(tdf))
        if i == 1
            Plots.plot(x, c, legend=false)
        else
            Plots.plot!(x, c, label = label[i])
        end
    end

    if TimeDataFrames.ncol(tdf) == 1 && isempty(title)
            title = names(tdf)[1]
    end
    Plots.plot!(label = label,
                title = title)
end
