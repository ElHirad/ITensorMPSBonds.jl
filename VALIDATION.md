# Validation — version 0.2.0, 2026-10-06

Tested on macOS arm64, Julia 1.11.5, and ITensors 0.9.32.

| ITensorMPS version | Tests passed | Execution |
| --- | ---: | --- |
| 0.3.45 | 135 / 135 | Complete test file with bounds checks in a separate compatibility environment |
| 0.4.2 | 135 / 135 | Complete test file with Julia's bounds checks enabled |

The independent reference uses only Julia arrays and LinearAlgebra SVDs. Tests
check reconstructed amplitudes against a separate TT-SVD implementation, then
reshape the complete small output at each cut to measure its actual Schmidt
rank. They cover Float32, Float64, ComplexF32, ComplexF64, integer arrays, views,
multidimensional arrays, ITensors, existing MPSs, tensor vectors, product-state
constructor forms, fixed-charge QN states, grouped sites, and dangling indices.

Failure tests include invalid dimensions, impossible neighboring-rank profiles,
rank-deficient product states, zero states, floating-point thresholds, and a GHZ
example where a later truncation lowers an earlier rank. Input preservation and
the final orthogonality center are checked too.

For version 0.1.1, a fresh temporary environment with `JULIA_LOAD_PATH=@:@stdlib` successfully
installed the package and its explicitly imported dependencies, ran the supplied
example, and confirmed that no Makie, CairoMakie, or Plots module was loaded.
Only shared dependency caches were reused. No upstream constructor methods are
overwritten; package precompilation and ambiguity checks passed.

The 0.2.0 test suite uses `exact_mps` throughout and loads both packages together,
checking that `MPS` still refers to the upstream constructor. The numerical
implementation and dependencies are unchanged from 0.1.1.

The included CI workflow covers Julia 1.10 and current stable Julia, and includes
Linux, macOS, and Windows jobs. The [previous revision passed remote CI](https://github.com/ElHirad/ITensorMPSBonds.jl/actions/runs/37520225814);
the 0.2.0 revision will trigger the same checks when pushed. Accelerators and arbitrary-precision backends were not tested;
their support depends on ITensors' conversion and SVD implementations.
