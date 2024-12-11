include("finite-elements.jl")

using .FiniteElements: gauss_quadrature_table
using .FiniteElements: generate_space

using .FiniteElements: build_LG, build_EQ
using .FiniteElements: phis_f, phi_derivs_f

using .FiniteElements: build_vec_2d
using .FiniteElements: build_mat_2d

using .FiniteElements: build_vec_mat_2d, build_vec_mat_2d_ref
using .FiniteElements: build_vec_mat_2d_iter, build_vec_mat_2d_iterref

include("bacarmo.jl")
using .Bacarmo: monta_K_quadrilatero, monta_F_quadrilatero

using BenchmarkTools: @benchmarkable, BenchmarkGroup
import BenchmarkTools as BT

const names = [
    :baseline, :ref, :iter, :iter_ref, :bacarmo,
]
const tests = [
    :build_vec, :build_mat, :build_both,
]
const tests_groups = Dict(
    :build_vec  => [
        (:baseline, build_vec_2d),
    ],
    :build_mat  => [
        (:baseline, build_mat_2d),
    ],
    :build_both => [
        (:baseline, build_vec_mat_2d),
        (:ref     , build_vec_mat_2d_ref),
        (:iter    , build_vec_mat_2d_iter),
        (:iter_ref, build_vec_mat_2d_iterref),
    ],
)
const tests_tags = [
    (:build_vec,  ["build_vec"]),
    (:build_mat,  ["build_mat"]),
    (:build_both, ["build_both"]),
]

const names_fullfuncs = (;
    baseline = build_vec_mat_2d,
    ref      = build_vec_mat_2d_ref,
    iter     = build_vec_mat_2d_iter,
    iter_ref = build_vec_mat_2d_iterref,
)
const names_tags = [
    (:baseline, ["no_ref", "no_iter"]),
    (:ref,      [   "ref", "no_iter"]),
    (:iter,     ["no_ref",    "iter"]),
    (:iter_ref, [   "ref",    "iter"]),
    (:bacarmo,  ["bacarmo"          ]),
]

suite = BenchmarkGroup()
for (test, test_tags) in tests_tags
    suite[test] = BenchmarkGroup(test_tags)
    for (name, name_tags) in names_tags
        suite[test][name] = BenchmarkGroup(name_tags)
    end
end

const f = (x...) -> foldl(+, x)
const alpha, beta = 1.0, 1.0

const min_max = 2:7

for i in min_max
    local pow = 1 << i
    local Ni = [ pow for _ in 1:2 ]
    local hi = 1 ./ Ni
    local N_e = foldl(*, Ni)

    local LG = build_LG(Ni)
    local m, EQ = build_EQ(Ni)
    local EQoLG = EQ[LG]

    local X, Y = generate_space(
        hi, Ni,
        noise = true,
        seed = 0,
    )

    for (name, func) in tests_groups[:build_vec]
        suite[:build_vec][name][pow] = @benchmarkable begin
            local gauss_n = 5
            local ws, ps = gauss_quadrature_table[gauss_n]

            local phis = phis_f(ps)
            local phi_derivs = phi_derivs_f(ps)

            F = $(func)(
                $f,
                $X, $Y, $N_e, $LG, $EQoLG, $m,
                phis, phi_derivs,
                ws, gauss_n,
            )
        end
    end
    suite[:build_vec][:bacarmo][pow] = @benchmarkable begin
        monta_F_quadrilatero($f, $X, $Y, $m, $EQ, $LG)
    end

    for (name, func) in tests_groups[:build_mat]
        suite[:build_mat][name][pow] = @benchmarkable begin
            local gauss_n = 5
            local ws, ps = gauss_quadrature_table[gauss_n]

            local phis = phis_f(ps)
            local phi_derivs = phi_derivs_f(ps)

            $(func)(
                $alpha, $beta,
                $X, $Y, $N_e, $LG, $EQoLG, $m,
                phis, phi_derivs,
                ws, gauss_n,
            )
        end
    end
    suite[:build_mat][:bacarmo][pow] = @benchmarkable begin
        monta_K_quadrilatero($alpha, $beta, $X, $Y, $m, $EQ, $LG)
    end

    for (name, func) in tests_groups[:build_both]
        suite[:build_both][name][pow] = @benchmarkable begin
            local gauss_n = 5
            local ws, ps = gauss_quadrature_table[gauss_n]

            local phis = phis_f(ps)
            local phi_derivs = phi_derivs_f(ps)

            $(func)(
                $f, $alpha, $beta,
                $X, $Y, $N_e, $LG, $EQoLG, $m,
                phis, phi_derivs,
                ws, gauss_n,
            )
        end
    end
    suite[:build_both][:bacarmo][pow] = @benchmarkable begin
        local K = monta_K_quadrilatero($alpha, $beta, $X, $Y, $m, $EQ, $LG)
        local F = monta_F_quadrilatero($f, $X, $Y, $m, $EQ, $LG)
        F, K
    end
end

BT.tune!(
    suite,
    verbose=true
)

bench = BT.run(
    suite,
    verbose=true,
)

BT.save("results.json", bench)
