# Security review — 2026-10-06

Scope: ITensorMPSBonds 0.1.1 source, examples, tests, repository history, CI
workflow, and the resolved dependency environment on macOS arm64 with Julia
1.11.5. This review is not a penetration test or a guarantee that the entire
Julia and native-library stack is free of vulnerabilities.

## Source and repository findings

No exploitable vulnerability was identified in this package's own source.

- Calculations use in-memory tensors and arrays. The package performs no network
  requests, file writes, subprocess execution, dynamic evaluation, deserialization,
  or input-controlled source inclusion. The three source includes are static.
- Dimension lengths, positive integer ranks, neighboring-rank feasibility, and
  numerical tolerances are validated. Feasibility products use arbitrary-size
  integers to avoid integer overflow in those checks.
- Inputs are copied where mutation would otherwise be possible. The final ranks
  are checked independently of the allocated link sizes. Numerical limitations
  are documented in the README rather than hidden by padding.
- Julia objects and callbacks are trusted executable input. Large calculations
  can exhaust resources; this package is not a server-side isolation boundary.
- Scanned all 12 original history blobs for common private-key, GitHub-token,
  AWS-key, and credential-assignment patterns; no matches found. This is a pattern
  scan, not an exhaustive secret-detection guarantee.
- CI uses immutable action commit pins, read-only `contents` permission, no
  persisted checkout credentials, no secrets, hosted runners, and a 30-minute
  timeout. It does not use `pull_request_target` or interpolate untrusted event
  text into executable commands.
- Public repository advisory APIs returned no published advisories for
  [ITensors](https://github.com/ITensor/ITensors.jl/security/advisories),
  [ITensorMPS](https://github.com/ITensor/ITensorMPS.jl/security/advisories),
  [actions/checkout](https://github.com/actions/checkout/security/advisories), or
  [setup-julia](https://github.com/julia-actions/setup-julia/security/advisories).
  Empty advisory lists do not establish that all upstream code is secure.

## Dependency comparison

The [version inventory](security-dependencies.csv) contains 106 packages. The
direct dependencies in the reviewed environment were ITensors 0.9.32 and
ITensorMPS 0.4.2. The source also passed tests with ITensorMPS 0.3.45.

The scan used 1,717 published advisories in the official
[Julia advisory database at dcb1dba](https://github.com/JuliaLang/SecurityAdvisories.jl/tree/dcb1dba21affb432702e3ec0e08d5753343c5a10).
That was the current repository revision when fetched on the review date; its
commit date is 2026-08-24, so this snapshot cannot establish coverage of every
later disclosure. The database also explicitly describes native-artifact
coverage as incomplete.

The scanner excludes withdrawn entries and compares affected ranges using Julia
version ordering, including JLL build suffixes. It stops on unsupported range
syntax. There were **77 advisory-range matches**, all in libraries bundled with
the old Julia 1.11.5 runtime. These are not 77 demonstrated exploits reachable
through this package. No matches were found for the separately installed package
dependencies in this database snapshot.

| Runtime component | Version | Range matches |
| --- | --- | ---: |
| LibCURL_jll | 8.6.0+0 | 37 |
| MbedTLS_jll | 2.28.6+0 | 17 |
| LibSSH2_jll | 1.11.0+1 | 11 |
| LibGit2_jll | 1.7.2+0 | 6 |
| nghttp2_jll | 1.59.0+0 | 3 |
| Zlib_jll | 1.2.13+1 | 2 |
| p7zip_jll | 17.4.0+2 | 1 |

The [complete matches](security-advisory-matches.csv) include upstream links and
descriptions. Several concern networking, Git operations, or archive extraction
rather than tensor arithmetic; exploitability through the package has not been
established. The package does not load a plotting or media-decoding stack.

**Required maintenance:** update Julia itself to the latest patch of a maintained
stable or LTS release, then resolve dependencies and repeat the scan for that
environment. `Pkg.update()` alone does not replace Julia's bundled libraries.
This review does not assert that an arbitrary Julia upgrade fixes every advisory.
The local runtime was not changed as part of publishing this package.

## Validation and reproduction

All 135 tests passed with each of ITensorMPS 0.3.45 and 0.4.2, including independent
dense-SVD comparisons, invalid inputs, input preservation, and exact-rank checks.
See [VALIDATION.md](../VALIDATION.md). Release 0.1.1 adds publication and security
documentation; the numerical implementation is unchanged from the tested 0.1.0.

To scan your actual environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
git clone https://github.com/JuliaLang/SecurityAdvisories.jl.git /path/to/advisories
julia --startup-file=no security/audit_dependencies.jl Manifest.toml /path/to/advisories report.toml
```

The included scanner parses advisory metadata as TOML data, does not execute it,
and does not send the manifest anywhere. It writes the requested report file;
choose a destination you control. Successful execution means that report
generation completed, **not** that the environment has no vulnerabilities.
Use current advisories for a new review; the pinned database commit above is for
reproducing this historical snapshot.
