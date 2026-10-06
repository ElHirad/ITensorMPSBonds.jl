function _from_tensor(a::ITensor, sites, chi;
        leftinds=nothing, orthocenter::Integer=length(sites),
        tags=["Link,l=$b" for b in 1:(length(sites) - 1)],
        rank_atol=0, rank_rtol=nothing, kwargs...)
    _check_keywords(kwargs)
    n = length(sites)
    n > 0 || throw(ArgumentError("Tensor decomposition requires at least one site."))
    1 <= orthocenter <= n || throw(ArgumentError("orthocenter must be a site number."))
    length(tags) == n - 1 || throw(DimensionMismatch("Provide one tag set per internal bond."))
    all(s -> hasinds(a, s), sites) || throw(ArgumentError("The input tensor must contain every site index."))
    isnothing(leftinds) || hasinds(a, leftinds) || throw(ArgumentError("leftinds are missing from the input tensor."))
    groups = [s isa Index ? [s] : collect(s) for s in sites]
    physical = reduce(vcat, groups; init=Index[])
    allunique(physical) || throw(ArgumentError("Site indices must not be repeated."))
    left = isnothing(leftinds) ? Index[] : leftinds isa Index ? [leftinds] : collect(leftinds)
    isempty(intersect(physical, left)) || throw(ArgumentError("leftinds must be distinct from site indices."))
    dims = [prod(big(dim(i)) for i in group; init=big(1)) for group in groups]
    dims[1] *= prod(big(dim(i)) for i in left; init=big(1))
    dims[end] *= prod(big(dim(i)) for i in uniqueinds(a, vcat(physical, left)); init=big(1))
    _feasible(chi, dims)
    tensors = ITensor[]
    remainder = a
    incoming = left
    for b in 1:(n - 1)
        rowinds = vcat(groups[b], incoming)
        u, remainder = _split(remainder, rowinds, chi[b], b;
            rank_atol, rank_rtol, tags=tags[b], kwargs...)
        incoming = [commonind(u, remainder)]
        push!(tensors, u)
    end
    push!(tensors, copy(remainder))
    return _verify_and_center(ITensorMPS.MPS(tensors; ortho_lims=n:n), chi;
        orthocenter, rank_atol, rank_rtol)
end

function _resize(input::ITensorMPS.MPS, chi;
        orthocenter::Integer=length(input), rank_atol=0, rank_rtol=nothing, kwargs...)
    _check_keywords(kwargs)
    m = deepcopy(input)
    n = length(m)
    n == 0 && return m
    1 <= orthocenter <= n || throw(ArgumentError("orthocenter must be a site number."))
    for b in 1:(n - 1)
        length(commoninds(m[b], m[b + 1])) == 1 || throw(ArgumentError(
            "Bond $b must have exactly one shared link index. Initialize the MPS with site indices first."))
    end
    _feasible(chi, [prod(big(dim(i)) for i in ITensorMPS.siteinds(m, b); init=big(1)) for b in 1:n])
    ITensorMPS.orthogonalize!(m, 1)
    for b in 1:(n - 1)
        old = ITensorMPS.linkind(m, b)
        left, remainder = _split(m[b], uniqueinds(m[b], m[b + 1]), chi[b], b;
            rank_atol, rank_rtol, tags=tags(old), kwargs...)
        m[b] = left
        m[b + 1] = remainder * m[b + 1]
    end
    m = ITensorMPS.MPS(collect(m); ortho_lims=n:n)
    return _verify_and_center(m, chi; orthocenter, rank_atol, rank_rtol)
end
