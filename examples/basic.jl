using ITensors
using ITensorMPS
using ITensorMPSBonds
using Random, LinearAlgebra

Random.seed!(14)
s = siteinds("Qubit", 6)
psi = randn(ComplexF64, 2^length(s))
psi /= norm(psi)
chi = [2, 3, 4, 3, 2]
m = exact_mps(psi, s, maxbond=chi)
@assert ITensorMPS.linkdims(m) == chi
println("Requested ranks: ", chi)
println("Returned dimensions: ", ITensorMPS.linkdims(m))
println("State error after truncation: ", norm(vec(Array(prod(m), s...)) - psi))

full = exact_mps(psi, s; maxbond=[2, 4, 8, 4, 2])
@assert vec(Array(prod(full), s...)) ≈ psi
println("Full-rank reconstruction agrees with the input.")

try
    exact_mps(s, "0"; maxbond=chi)
    error("A product state must not acquire artificial nonzero ranks.")
catch e
    e isa ArgumentError || rethrow()
    println("Impossible product-state request correctly rejected: ", sprint(showerror, e))
end
