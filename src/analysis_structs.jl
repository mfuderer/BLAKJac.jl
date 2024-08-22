"""
    struct TrajectoryElement{T<:Real}

Representing a profile of a trajectory with certain `ky` and `kz` indices.
"""
struct TrajectoryElement{T<:Real}
    ky::T
    kz::T
end

"""
    TrajectorySet(ky::Vector{Float64}, kz::Vector{Float64})

Converts two arrays ky and kz into a set of (single-element) 'trajectories'

The arrays have to be of the same length: the number of profiles acquired during the sequence. If 2D sequences are considered, kz 
can be initialized by kz=ones(nTR). The values can be provided 'signed' (e.g. ky ranging from -112 to +111) or 'unsigned' 
(e.g. ranging from 1 to 224). In the latter case, the range of each array is first translated to a symmetric range around 0.
"""
function TrajectorySet(ky::Vector{T}, kz::Vector{T}) where T<:Real

    @assert length(ky) == length(kz)
    # some tricky logic follows: if any of the provided vectors shows negative values, it is assumed that the k-space center 
    # corresponds to the zero-value of input; if the input is nonnegative, it is assumed that it is symmetric around k=0
    ky = _to_signed_k_range(ky)
    kz = _to_signed_k_range(kz)

    trajectorySet = [[TrajectoryElement(ky[i], kz[i])] for i in 1:length(ky)]

    return trajectorySet
end

function _to_signed_k_range(k::Vector{<:Real})
    # Don't do anything if already 'signed'
    (minimum(k) < 0) && return k
    # Calcalate the 'unsigned' center value
    ck = floor(maximum(k) / 2) + 1.0
    # Translate the k values to 'signed' values with the kspace center at 0
    return k .- ck
end

"""
    nTR(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})

Extracts the number of profiles from the trajectorySet.
"""
function nTR(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    return length(trajectorySet)
end

"""
    nky(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})

Extracts the number of ky values per k-space from the trajectorySet. Assumes that the outer ky values are part of trajectorySet.
"""
function nky(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    ky = [s[i].ky for s in trajectorySet for i in 1:length(s)]
    return round(Int64, maximum(ky) - minimum(ky) + 1.0)
end

"""
    nkz(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})

Extracts the number of kz values per k-space from the trajectorySet. Assumes that the outer kz values are part of trajectorySet.
"""
function nkz(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    kz = [s[i].kz for s in trajectorySet for i in 1:length(s)]
    return round(Int64, maximum(kz) - minimum(kz) + 1.0)
end

