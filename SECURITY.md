# Security policy

Use a maintained Julia release with current security fixes and current compatible
package dependencies. Julia version compatibility does not imply that an older
runtime is still secure or maintained.

## Intended use

ITensorMPSBonds operates on in-memory scientific data and returns an MPS. Its
source has no network access, file writing, subprocess execution, dynamic code
evaluation, or deserialization. It does not modify upstream constructor methods.

Julia functions, custom array implementations, and custom numeric types run with
the caller's permissions. This library is not a sandbox for untrusted Julia code.
Large inputs and SVD computations can consume substantial memory and CPU. Services
accepting external inputs must impose their own size limits and process isolation.
The library's dimension validation is not a resource quota.

The [security review](docs/SECURITY_REVIEW.md) records the scope, tested dependency
versions, checks performed, remaining advisory matches, and coverage limitations.
No review can guarantee the absence of undiscovered vulnerabilities.

## Reporting

Report exploitable issues privately to the maintainer or the affected upstream
project using its security policy. If a private contact is unavailable, the
[Julia security team](https://security.julialang.org/about/) accepts ecosystem
reports at `security@julialang.org`. Do not post credentials, private data, or
working exploit details in a public issue. Include versions, a minimal
reproduction, and the affected operation, with sensitive information removed.
