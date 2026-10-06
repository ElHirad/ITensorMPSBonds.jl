# ITensorMPSBonds.jl

Construct a standard `ITensorMPS.MPS` with a specified **nonzero Schmidt rank at
every internal bond**. Uses ITensors SVDs, accepts arrays and ITensors, and supports
the ordinary ITensorMPS product-state constructor forms. No plotting dependencies.

```julia
using ITensors
using ITensorMPS
using ITensorMPSBonds

s = siteinds("Qubit", 6)
psi = randn(ComplexF64, 64)
chi = [2, 3, 4, 3, 2]  # N sites have N-1 internal bonds

m = exact_mps(psi, s, maxbond=chi)
@assert m isa ITensorMPS.MPS
@assert ITensorMPS.linkdims(m) == chi
```

The keyword spelling `exact_mps(psi, s; maxbond=chi)` is equivalent. The output works
with ITensorMPS operations and the other ITensor packages.

## Installation

Requires Julia 1.10 or later and a maintained Julia runtime is recommended.
This standalone package is not registered in General. Install from the public
GitHub repository once in the Julia environment used by your notebook or project:

```julia
using Pkg
Pkg.add(url="https://github.com/ElHirad/ITensorMPSBonds.jl", rev="v0.2.0")
Pkg.add(["ITensors", "ITensorMPS"])  # For the explicit imports in your own code
```

The release tag selects a fixed package version. Omit `rev` to track the main
branch instead. If you previously installed a local development copy, the same
`Pkg.add(url=..., rev=...)` call switches it to the GitHub release.

Julia resolves the dependencies automatically; the second line makes them direct
dependencies of your environment so you can import them yourself. The first import may compile dependencies;
normal scripts need only the imports shown above. This package does not load
Makie, CairoMakie, or either of the previous plotting/drawing packages.

For local development or an extracted ZIP, use
`Pkg.develop(path="/absolute/path/to/ITensorMPSBonds")` and keep that folder in place.

## Exact ranks, not zero padding

`maxbond` is a vector of exactly N-1 positive integers. Despite the keyword name,
the requested dimensions are **exact**, not upper bounds. There is no zero padding
or added random noise. A successful call checks both link sizes and nonzero
Schmidt values; impossible requests throw an `ArgumentError`.

- Reducing a rank generally changes the represented state. The method performs
  sequential SVD truncation; it is not a globally optimal variational fit.
- The default `cutoff=0.0` keeps up to the requested number of values at each cut.
  The state is not automatically normalized. Use `ITensorMPS.normalize!(m)` if needed.
- A product state has rank one at every cut, so requesting a larger rank throws.
- Rank profiles must respect physical and neighboring-rank bounds. For example,
  `[2, 4, 1]` is impossible for four qubits, even though 4 is the unconstrained
  maximum at the middle cut.
- Truncating a later cut can lower an earlier rank. A final canonical SVD sweep
  checks the completed state and rejects it if that happens. Rejection means
  this construction did not achieve the ranks; it is not proof that no other,
  differently chosen approximation could achieve them.

In floating-point arithmetic, “nonzero” needs a numerical threshold. A singular
value counts only if it exceeds
`max(rank_atol, rank_rtol * largest_singular_value)`. Defaults are `rank_atol=0`
and `rank_rtol=max(local SVD matrix dimensions)*eps(real scalar type)` at each
check. These checks include the final canonical sweep. Override deliberately
when your problem has meaningful tiny coefficients:

```julia
m = exact_mps(psi, s; maxbond=chi, rank_rtol=1e-12, rank_atol=0.0)
```

Setting both tolerances to zero uses strict positivity and can count floating-point
roundoff as rank. An exactly zero state has no nonzero Schmidt values and is
rejected for any chain with internal bonds.

## Inputs and compatibility

| Input | Example |
| --- | --- |
| Numeric vector, matrix, multidimensional array, or array view | `exact_mps(psi, s; maxbond=chi)` |
| Dense or QN-conserving ITensor | `exact_mps(A, s; maxbond=chi)` |
| Existing MPS | `exact_mps(m; maxbond=chi)` or `exact_mps(m, s; maxbond=chi)` |
| Vector of site ITensors | `exact_mps(collect(m); maxbond=chi)` |
| Product-state string or integer labels | `exact_mps(s, "Up"; maxbond=ones(Int, length(s)-1))` |
| Product-state label vector or function | `exact_mps(s, n -> isodd(n) ? "Up" : "Dn"; maxbond=ones(Int, length(s)-1))` |
| Explicit product-state scalar type | `exact_mps(ComplexF32, s, "Up"; maxbond=ones(Int, length(s)-1))` |
| Index-value pairs, with optional scalar type | `exact_mps([i => 1 for i in s]; maxbond=ones(Int, length(s)-1))` |

Array conversion and ordering follow `ITensorMPS.MPS`: Julia column-major order,
with the first supplied site varying fastest. Real/complex single/double precision
are supported without forced conversion to Float64. Other scalar types and array
backends retain the capabilities and restrictions of ITensors' conversion and SVD
implementations; arbitrary precision and accelerators are not independently tested.

QN inputs retain symmetry and flux. The integer dimensions specify total ranks,
not how many states to keep in each charge sector; ITensors selects the retained
singular values. A sector-restricted state can have fewer available nonzero values
than the dense physical bound and is then rejected.

An ITensor may also use grouped site indices or dangling boundary indices, following
the upstream convention. `leftinds`, per-bond `tags`, and `orthocenter` are supported
for that decomposition. Existing MPS inputs can choose `orthocenter` too; the default
is the last site. Inputs are not mutated.

SVD options such as `svd_alg`, `cutoff`, `use_relative_cutoff`, and `min_blockdim`
are forwarded for array, ITensor, and existing-MPS inputs. Options that prevent
the requested exact ranks cause rejection. Conflicting controls such as `maxdim`,
`mindim`, and `linkdims` are rejected when using `maxbond`.

Without `maxbond`, all arguments and ordinary constructor keywords are forwarded
to `ITensorMPS.MPS`, including empty/uninitialized constructors. With `maxbond`,
the input must describe an initialized state; an empty tensor allocation cannot
satisfy positive nonzero ranks. A single-site state uses `maxbond=Int[]`.

## Constructor name and upgrading

The package exports only `exact_mps`. You can load `ITensorMPS` and
`ITensorMPSBonds` together: `MPS(...)` remains the ordinary upstream constructor,
while `exact_mps(...; maxbond=chi)` constructs the requested exact ranks.
The returned **type** is always `ITensorMPS.MPS`, so existing MPS operations
and type annotations continue to work. No upstream methods are overwritten.

Version 0.2.0 removes the package's previous `MPS` name. When upgrading from
0.1.x, replace `using ITensorMPSBonds: MPS` with `using ITensorMPSBonds`, and
replace calls to this package's constructor with `exact_mps(...)`.
Calls to the ordinary `ITensorMPS.MPS(...)` constructor do not change.
If you already used `exact_mps`, no call changes are needed.
Restart an existing notebook kernel after upgrading to clear old imports.

## Examples and verification

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia --project=. examples/basic.jl
```

Tests compare the reconstruction against an independent dense TT-SVD, check
actual ranks by reshaping the complete output at every cut, and cover scalar
types, product-state forms, QN conservation, input preservation, tolerance
boundaries, and impossible requests. Dense reconstruction is used only in small
tests, never to resize an existing MPS in the package.
The final rank verification adds a canonical SVD sweep. See
[VALIDATION.md](VALIDATION.md) for the tested versions and results.

## Security

See [SECURITY.md](SECURITY.md) and the [dated security review](docs/SECURITY_REVIEW.md).
The package has no network, file-writing, shell-execution, or deserialization
feature. Use it with trusted Julia code and keep Julia itself up to date: package
updates alone do not replace vulnerable libraries bundled with an older runtime.

## License

MIT. The implementation builds on the public
[ITensors](https://github.com/ITensor/ITensors.jl) and
[ITensorMPS](https://github.com/ITensor/ITensorMPS.jl) APIs.
