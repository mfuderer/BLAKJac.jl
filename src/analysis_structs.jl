struct TrajectoryElement
    ky::Float64
    kz::Float64
end

"""
    TrajectorySet(ky::Vector{Float64}, kz::Vector{Float64})

Converts two arrays ky and kz into a set of (single-element) 'trajectories'

The arrays have to be of the same length: the number of profiles acquired during the sequence. If 2D sequences are considered, kz 
can be initialized by kz=ones(nTR). The values can be provided 'signed' (e.g. ky ranging from -112 to +111) or 'unsigned' 
(e.g. ranging from 1 to 224). In the latter case, the range of each array is first translated to a symmetric range around 0.
"""
function TrajectorySet(ky::Vector{Float64}, kz::Vector{Float64})
    @assert length(ky) == length(kz)
    # some tricky logic follows: if any of the provided vectors shows negative values, it is assumed that the k-space center 
    # corresponds to the zero-value of input; if the input is nonnegative, it is assumed that it is symmetric around k=0
    already_signed = (minimum(ky) < 0 || minimum(kz) < 0)
    cky = already_signed ? 0.0 : floor(maximum(ky) / 2) + 1.0
    ckz = already_signed ? 0.0 : floor(maximum(kz) / 2) + 1.0
    trajectorySet = Vector{Vector{TrajectoryElement}}(undef, length(ky))
    for i in 1:length(ky)
        tr = TrajectoryElement(ky[i] - cky, kz[i] - ckz)
        trajectorySet[i] = Vector{TrajectoryElement}([tr])
    end
    return trajectorySet
end

# "conditional conjugation"
function CondConj(cond::Bool, f::ComplexF64)
    cf = cond ? conj(f) : f
    return cf
end
function CondConj(cond::Bool, f::Vector{ComplexF64})
    cf = cond ? conj.(f) : f
    return cf
end