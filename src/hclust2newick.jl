#=
    transliterated from  
        https://rdrr.io/bioc/ctc/src/R/hc2Newick.R
    2021 spring  TEC
=#

function hclust2newick(hclu::Hclust{Float64}; labels::Any = nothing, flat::Bool= true)

    if isnothing(labels)
        labels = 1:length(hclu.order)
    # else
    #   labels = hclu.labels   Hclust.labels flagged for deprecation
    end

    function putparenthesis(i)
        ## recursive function
        ## i: index of row in hclu.merge
        j = hclu.merge[i, 1]
        k = hclu.merge[i, 2]
        dist = 0.0

        if j < 0
            left = labels[-j]
            if k > 0
                dist = hclu.height[i] - hclu.height[k]
            else
                dist = hclu.height[i]
            end
        else
            left = putparenthesis(j)
        end

        if k < 0
            right = labels[-k]
            if j > 0
                dist = hclu.height[i] - hclu.height[j]
            else
                dist = hclu.height[i]
            end
        else
            right = putparenthesis(k)
        end

        if flat
            return join( # R paste()
                ["(", left, ":", dist/2.0, ",", right, ":", dist/2.0, ")"], "")  
        else  # I do not know this use case. so dict of dicts it is for now
            return Dict(:left=>left, :right=>right, :dist=>dist)
        end
    end

    newick = putparenthesis(first(size(hclu.merge)))  # R nrows
    if flat
        newick = join([newick, ";"], "")  # R paste()
    end
    return newick
end


function write_newick(io::IO, hclu::Hclust{Float64}; labels::Any = nothing)
    write(io, hclust2newick(hclu, labels = labels, flat = true))
end

