using Test, Random, LinearAlgebra
using ITensors
import ITensorMPS
using ITensorMPS: siteinds
import ITensorMPSBonds
using ITensorMPSBonds: MPS, exact_mps

Random.seed!(127)
dense(m, s) = vec(Array(prod(m), s...))
function dense_ranks(m, s; rtol=1e-10)
    a = dense(m, s)
    return [count(>(rtol * maximum(sv)), sv) for b in 1:(length(s) - 1)
        for sv in (svdvals(reshape(a, prod(dim, s[1:b]), :)),)]
end

# Independent array-only TT-SVD reference. Never calls ITensors or package internals.
function reference_tt(a, dims, chi)
    cores = Array[]
    remainder = reshape(copy(a), 1, :)
    r = 1
    for b in eachindex(chi)
        f = LinearAlgebra.svd(reshape(remainder, r * dims[b], :))
        push!(cores, reshape(f.U[:, 1:chi[b]], r, dims[b], chi[b]))
        remainder = Diagonal(f.S[1:chi[b]]) * f.Vt[1:chi[b], :]
        r = chi[b]
    end
    push!(cores, reshape(remainder, r, dims[end], 1))
    out = reshape(cores[1], dims[1], chi[1])
    for b in 2:length(dims)
        out = reshape(out * reshape(cores[b], size(cores[b], 1), :), :, size(cores[b], 3))
    end
    return vec(out)
end

@testset "Per-cut TT-SVD and dense reconstruction" begin
    for T in (Float32, Float64, ComplexF32, ComplexF64)
        s = siteinds("Qubit", 6)
        a = randn(T, 64)
        chi = [2, 3, 4, 3, 2]
        original = copy(a)
        m = MPS(a, s, maxbond=chi)
        @test m isa ITensorMPS.MPS
        @test ITensorMPS.linkdims(m) == chi
        @test eltype(m[1]) == T
        @test a == original
        tol = T in (Float32, ComplexF32) ? 1e-5 : 1e-12
        @test dense(m, s) ≈ reference_tt(a, fill(2, 6), chi) rtol=tol
        @test dense_ranks(m, s; rtol=10tol) == chi
        full = MPS(a, s; maxbond=[2, 4, 8, 4, 2])
        @test dense(full, s) ≈ a rtol=tol
        @test norm(full) ≈ norm(a) rtol=tol
        @test ITensorMPS.isortho(m)
        @test ITensorMPS.orthocenter(m) == 6
    end
end

@testset "Array representations and ITensor inputs" begin
    s = siteinds("Qubit", 4)
    a = randn(16)
    inputs = (a, reshape(a, 4, 4), reshape(a, 2, 2, 2, 2), view(vcat(a, a), 1:16),
        PermutedDimsArray(reshape(a, 4, 4), (2, 1)), round.(Int, 100a))
    for data in inputs
        native = ITensorMPS.MPS(data, s; cutoff=0)
        m = MPS(data, s; maxbond=[2, 4, 2])
        @test dense(m, s) ≈ dense(native, s)
    end
    tensor = ITensor(a, s...)
    original = copy(tensor)
    m = MPS(tensor, s; maxbond=[2, 3, 2], orthocenter=2, tags=["MyLink,b=$b" for b in 1:3])
    @test dense_ranks(m, s) == [2, 3, 2]
    @test ITensorMPS.orthocenter(m) == 2
    @test tensor == original
    @test all(b -> hastags(ITensorMPS.linkind(m, b), "MyLink,b=$b"), 1:3)
    @test exact_mps === MPS
end

@testset "Existing MPS and tensor vectors" begin
    s = siteinds("Qubit", 6)
    a = randn(ComplexF64, 64)
    original = ITensorMPS.MPS(a, s; cutoff=0)
    saved = deepcopy(original)
    for input in (original, collect(original))
        m = MPS(input; maxbond=[2, 3, 4, 3, 2])
        @test ITensorMPS.linkdims(m) == [2, 3, 4, 3, 2]
        @test dense_ranks(m, s) == [2, 3, 4, 3, 2]
        @test dense(m, s) ≈ reference_tt(a, fill(2, 6), [2, 3, 4, 3, 2])
    end
    m = MPS(original, s; maxbond=[2, 4, 8, 4, 2])
    @test dense(m, s) ≈ a
    m[1] *= 3
    @test all(b -> original[b] == saved[b], 1:6)
    @test ITensorMPS.orthocenter(original) == ITensorMPS.orthocenter(saved)
    @test_throws ArgumentError MPS(original, reverse(s); maxbond=[2, 4, 8, 4, 2])
end

@testset "Product state constructor forms" begin
    s = siteinds("S=1/2", 4)
    for state in ("Up", 1, ["Up", "Dn", "Up", "Dn"], [1, 2, 1, 2], n -> isodd(n) ? "Up" : "Dn")
        for args in ((s, state), (ComplexF32, s, state))
            m = MPS(args...; maxbond=[1, 1, 1])
            native = ITensorMPS.MPS(args...)
            @test dense(m, s) ≈ dense(native, s)
            @test ITensorMPS.linkdims(m) == [1, 1, 1]
        end
        @test_throws ArgumentError MPS(s, state; maxbond=[2, 2, 2])
    end
    for args in (([s[b] => 1 for b in 1:4],), (ComplexF64, [s[b] => 1 for b in 1:4]))
        @test ITensorMPS.linkdims(MPS(args...; maxbond=[1, 1, 1])) == [1, 1, 1]
    end
end

@testset "Impossible nonzero ranks and tolerance" begin
    s = siteinds("Qubit", 3)
    product = zeros(8); product[1] = 1
    ghz = zeros(8); ghz[1] = sqrt(0.8); ghz[end] = sqrt(0.2)
    @test_throws ArgumentError MPS(product, s; maxbond=[2, 2])
    @test_throws ArgumentError MPS(zeros(8), s; maxbond=[1, 1])
    @test_throws ArgumentError MPS(randn(8), s; maxbond=[3, 2])
    @test_throws ArgumentError MPS(randn(16), siteinds("Qubit", 4); maxbond=[2, 4, 1])
    # Both requested ranks pass locally, but the last truncation makes a product state.
    @test_throws ArgumentError MPS(ghz, s; maxbond=[2, 1])
    @test dense_ranks(MPS(ghz, s; maxbond=[2, 2]), s) == [2, 2]
    s2 = siteinds("Qubit", 2)
    tiny = [1.0, 0.0, 0.0, 1e-10]
    @test_throws ArgumentError MPS(tiny, s2; maxbond=[2], rank_rtol=1e-8)
    @test ITensorMPS.linkdims(MPS(tiny, s2; maxbond=[2], rank_rtol=1e-12)) == [2]
    @test_throws ArgumentError MPS(tiny, s2; maxbond=[2], rank_atol=1e-9)
end

@testset "QN conservation and rank checks" begin
    s = siteinds("S=1/2", 4; conserve_qns=true)
    # A fixed-charge state with maximal ranks [2,4,2].
    a = ITensor(ComplexF64, s...)
    for bits in Iterators.product(fill(1:2, 4)...)
        count(==(1), bits) == 2 || continue
        a[(s[b] => bits[b] for b in 1:4)...] = randn(ComplexF64)
    end
    before = copy(a)
    for chi in ([2, 4, 2], [2, 3, 2])
        m = MPS(a, s; maxbond=chi)
        @test hasqns(m)
        @test flux(m) == flux(a)
        @test ITensorMPS.linkdims(m) == chi
        @test dense_ranks(m, s) == chi
        if chi == [2, 4, 2]
            @test prod(m) ≈ a
        end
    end
    @test a == before
    product = MPS(s, ["Up", "Dn", "Up", "Dn"]; maxbond=[1, 1, 1])
    @test flux(product) == QN("Sz", 0)
    @test_throws ArgumentError MPS(product; maxbond=[2, 2, 2])
end

@testset "Grouped sites and dangling boundary indices" begin
    s = siteinds("Qubit", 4)
    a = random_itensor(s...)
    m = MPS(a, [(s[1], s[2]), (s[3], s[4])]; maxbond=[3])
    @test ITensorMPS.linkdims(m) == [3]
    @test Array(prod(m), s...) ≈ reshape(reference_tt(vec(Array(a, s...)), [4, 4], [3]), 2, 2, 2, 2)
    l, r = Index(2, "Left"), Index(2, "Right")
    a = random_itensor(l, s[1], s[2], r)
    m = MPS(a, s[1:2]; leftinds=l, maxbond=[4])
    @test prod(m) ≈ a
    @test hasind(m[1], l)
    @test hasind(m[2], r)
end

@testset "Invalid input and forwarding" begin
    s = siteinds("Qubit", 3)
    for chi in ([1], [1, 1, 1], Int[])
        @test_throws DimensionMismatch MPS(randn(8), s; maxbond=chi)
    end
    for chi in ([0, 1], [-1, 1], [1.0, 1.0], [true, true], [1, big(typemax(Int)) + 1])
        @test_throws ArgumentError MPS(randn(8), s; maxbond=chi)
    end
    @test_throws ArgumentError MPS(randn(8), s; maxbond=2)
    @test_throws ArgumentError MPS(randn(8), s; maxbond=[1, 1], maxdim=1)
    @test_throws ArgumentError MPS(randn(8), s; maxbond=[1, 1], cutoff=-1)
    @test_throws ArgumentError MPS(randn(8), s; maxbond=[1, 1], rank_rtol=-1)
    @test_throws ArgumentError MPS(randn(8), s; maxbond=[1, 1], rank_atol=NaN)
    @test_throws ArgumentError MPS(randn(8), s; maxbond=[1, 1], orthocenter=0)
    @test_throws ArgumentError MPS(random_itensor(s...), [s[1], s[1], s[3]]; maxbond=[1, 1])
    one = MPS([2.0, 3.0], s[1:1]; maxbond=Int[])
    @test dense(one, s[1:1]) ≈ [2.0, 3.0]
    @test isempty(MPS())
    @test length(MPS(3)) == 3
    @test ITensorMPS.linkdims(MPS(s, "0")) == [1, 1]
    @test MPS(randn(8), s; maxdim=2) isa ITensorMPS.MPS
    @test isempty(Test.detect_ambiguities(ITensorMPSBonds; recursive=true))
end
