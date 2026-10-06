function _bond_dimensions(maxbond, n::Integer)
    maxbond isa AbstractVector || throw(ArgumentError("maxbond must be a vector of positive integers."))
    length(maxbond) == max(0, n - 1) || throw(DimensionMismatch(
        "An MPS with $n sites needs $(max(0, n - 1)) internal bond dimensions."))
    all(x -> x isa Integer && !(x isa Bool) && 0 < x <= typemax(Int), maxbond) ||
        throw(ArgumentError("Every maxbond entry must be a positive integer representable as Int."))
    return Int[x for x in maxbond]
end

function _check_keywords(kwargs)
    for key in (:maxdim, :mindim, :linkdims, :ortho, :which_decomp, :eigen_perturbation)
        haskey(kwargs, key) && throw(ArgumentError(
            "$key conflicts with exact bond dimensions; use maxbond and the SVD options instead."))
    end
end

function _check_tolerances(rank_atol, rank_rtol)
    rank_atol isa Real && isfinite(rank_atol) && rank_atol >= 0 ||
        throw(ArgumentError("rank_atol must be finite and nonnegative."))
    isnothing(rank_rtol) || (rank_rtol isa Real && isfinite(rank_rtol) && rank_rtol >= 0) ||
        throw(ArgumentError("rank_rtol must be finite and nonnegative, or nothing."))
end

function _feasible(chi, physicaldims)
    n = length(physicaldims)
    for b in eachindex(chi)
        left = b == 1 ? 1 : chi[b - 1]
        right = b == n - 1 ? 1 : chi[b + 1]
        bound = min(big(left) * physicaldims[b], big(right) * physicaldims[b + 1])
        chi[b] <= bound || throw(ArgumentError(
            "Requested rank $(chi[b]) at bond $b is impossible with these site dimensions and neighboring ranks (upper bound $bound)."))
    end
end

function _require_rank(s::ITensor, target, b, matrixsize; rank_atol, rank_rtol)
    # Transfer only the diagonal for rank checks; avoid scalar indexing of
    # accelerator tensors and avoid materializing a dense chi-by-chi matrix.
    values = abs.(vec(Array(ITensors.diag(s))))
    all(isfinite, values) || throw(ArgumentError("Nonfinite singular values at bond $b."))
    rtol = isnothing(rank_rtol) ? max(matrixsize...) * eps(typeof(float(real(one(eltype(s)))))) : rank_rtol
    threshold = max(rank_atol, rtol * maximum(values; init=zero(eltype(values))))
    numerical_rank = count(>(threshold), values)
    numerical_rank >= target || throw(ArgumentError(
        "Bond $b requires $target nonzero Schmidt values, but only $numerical_rank exceed the numerical-rank threshold $threshold. No zero padding is performed. Reduce maxbond or, for intentionally tiny values, adjust rank_rtol/rank_atol."))
    length(values) == target || throw(ArgumentError(
        "SVD returned $(length(values)) values at bond $b instead of $target; check conflicting SVD options."))
    return nothing
end

function _split(a, rowinds, target, b; rank_atol, rank_rtol, tags,
        cutoff=0.0, svd_alg=nothing, kwargs...)
    isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError("cutoff must be finite and nonnegative."))
    sizes = (dim(commoninds(a, rowinds)), dim(uniqueinds(a, rowinds)))
    u, s, v = svd(a, rowinds; maxdim=target, mindim=1, cutoff,
        lefttags=tags, righttags=tags, alg=svd_alg, kwargs...)
    _require_rank(s, target, b, sizes; rank_atol, rank_rtol)
    return u, s * v
end

function _verify_and_center(m, chi; orthocenter, rank_atol, rank_rtol)
    n = length(m)
    n == 0 && return m
    1 <= orthocenter <= n || throw(ArgumentError("orthocenter must be a site number."))
    ITensorMPS.orthogonalize!(m, 1)
    # Later truncations can reduce earlier ranks. Check the completed state's
    # Schmidt values in a canonical sweep, rather than just inspecting link sizes.
    for b in eachindex(chi)
        left, remainder = _split(m[b], uniqueinds(m[b], m[b + 1]), chi[b], b;
            rank_atol, rank_rtol, tags=tags(ITensorMPS.linkind(m, b)))
        m[b] = left
        m[b + 1] = remainder * m[b + 1]
    end
    m = ITensorMPS.MPS(collect(m); ortho_lims=n:n)
    ITensorMPS.orthogonalize!(m, orthocenter)
    ITensorMPS.linkdims(m) == chi || error("Internal error: exact bond dimensions were not achieved.")
    return m
end
