include("finite-elements.jl")

using .FiniteElements: gauss_quadrature_table
using .FiniteElements: generate_space

using .FiniteElements: build_LG, build_EQ
using .FiniteElements: phis_f, phi_derivs_f

using .FiniteElements: build_vec_mat_2d, build_vec_mat_2d_ref
using .FiniteElements: build_vec_mat_2d_iter, build_vec_mat_2d_iterref

using BenchmarkTools: @benchmarkable, BenchmarkGroup
import BenchmarkTools as BT

const names_funcs = [
    (:baseline, build_vec_mat_2d),
    (:ref, build_vec_mat_2d_ref),
    (:iter, build_vec_mat_2d_iter),
    (:iter_ref, build_vec_mat_2d_iterref),
]
const tags = [
    ["no_ref", "no_iter"],
    [   "ref", "no_iter"],
    ["no_ref",    "iter"],
    [   "ref",    "iter"],
]

suite = BenchmarkGroup()
for i in 1:length(names_funcs)
    suite[names_funcs[i][1]] = BenchmarkGroup(tags[i])
end

const gauss_n = 5
const ws, ps = gauss_quadrature_table[gauss_n]

const phis = phis_f(ps)
const phi_derivs = phi_derivs_f(ps)

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

    for (name, func) in names_funcs
        suite[name][pow] = @benchmarkable $(func)(
            $f, $alpha, $beta,
            $X, $Y, $N_e, $LG, $EQoLG, $m,
            $phis, $phi_derivs,
            $ws, $gauss_n,
        )
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
