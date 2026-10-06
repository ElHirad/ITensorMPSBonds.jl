"""
    exact_mps(psi, sites; maxbond=[chi_1, ..., chi_Nminus1], kwargs...)
    exact_mps(args...; maxbond, kwargs...)

Create a standard `ITensorMPS.MPS` with exactly the requested internal link
dimensions. For N sites, supply N-1 positive integers. Tensor and array inputs
are decomposed with ITensors SVDs, limiting each cut individually. Existing MPSs
are copied and compressed by a left-to-right sweep. Every requested dimension
must correspond to nonzero Schmidt values. Impossible
requests are rejected; no zero padding is performed. The final ranks are verified
after all truncations, using configurable numerical-rank tolerances.

The default `cutoff=0.0` retains values up to the requested dimensions. Smaller
dimensions generally approximate the input state; no normalization is applied.
QN tensors keep their charges and block structure. Integers alone do not specify
a desired sector distribution. Rank uses `max(rank_atol, rank_rtol * sigma_max)`;
by default `rank_atol=0` and `rank_rtol=max(local matrix dimensions)*eps(scalar type)`.

Ordinary ITensorMPS constructor forms are forwarded, including product-state
labels, functions, explicit scalar types, Index-value pairs, and tensor vectors.
Without `maxbond`, arguments are passed straight to the original constructor.
`exact_mps(existing_mps, sites; maxbond)` additionally checks the site ordering.

Load this constructor with `using ITensorMPSBonds`. It can be used alongside
`using ITensorMPS`: the upstream `MPS` constructor keeps its original name.
The returned type is `ITensorMPS.MPS`; use that type for annotations.
No upstream methods are overwritten.
"""
function exact_mps(args...; maxbond=nothing, rank_atol=0, rank_rtol=nothing, kwargs...)
    isnothing(maxbond) && return ITensorMPS.MPS(args...; kwargs...)
    _check_keywords(kwargs)
    _check_tolerances(rank_atol, rank_rtol)
    if length(args) == 2 && args[1] isa ITensor
        a, sites = args
        chi = _bond_dimensions(maxbond, length(sites))
        return _from_tensor(a, sites, chi; rank_atol, rank_rtol, kwargs...)
    elseif length(args) == 2 && args[1] isa AbstractArray && !(args[1] isa AbstractVector{<:Index})
        a, sites = args
        chi = _bond_dimensions(maxbond, length(sites))
        return _from_tensor(itensor(a, sites...), sites, chi; rank_atol, rank_rtol, kwargs...)
    elseif length(args) in (1, 2) && args[1] isa ITensorMPS.MPS
        input = args[1]
        chi = _bond_dimensions(maxbond, length(input))
        if length(args) == 2
            sites = args[2]
            length(sites) == length(input) || throw(DimensionMismatch("Site count does not match the MPS."))
            all(b -> Set(sites[b] isa Index ? [sites[b]] : sites[b]) == Set(ITensorMPS.siteinds(input, b)), eachindex(sites)) ||
                throw(ArgumentError("Site indices must match the existing MPS in the supplied order."))
        end
        return _resize(input, chi; rank_atol, rank_rtol, kwargs...)
    end
    # Constructor-specific keywords belong to upstream; resizing defaults are
    # separate so product-state constructors need not implement SVD keywords.
    input = ITensorMPS.MPS(args...; kwargs...)
    chi = _bond_dimensions(maxbond, length(input))
    return _resize(input, chi; rank_atol, rank_rtol)
end
