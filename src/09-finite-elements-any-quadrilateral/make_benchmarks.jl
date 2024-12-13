include("finite-elements.jl")

using .FiniteElements: gauss_quadrature_table
using .FiniteElements: generate_space

using .FiniteElements: build_LG, build_EQ
using .FiniteElements: phis_f, phi_derivs_f
using .FiniteElements: x2xis_f, dx2xis_f

using .FiniteElements: build_small_vec_2d, build_small_vec_2d_ref!
using .FiniteElements: build_small_mat_2d, build_small_mat_2d_ref!

using .FiniteElements: build_vec_2d, build_vec_2d_ref
using .FiniteElements: build_vec_2d_iter, build_vec_2d_iterref

using .FiniteElements: build_mat_2d, build_mat_2d_ref
using .FiniteElements: build_mat_2d_iter, build_mat_2d_iterref

using .FiniteElements: build_vec_mat_2d, build_vec_mat_2d_ref
using .FiniteElements: build_vec_mat_2d_iter, build_vec_mat_2d_iterref

include("bacarmo.jl")
using .Bacarmo: monta_Kᵉ_quadrilatero!, monta_Fᵉ_quadrilatero!
using .Bacarmo: monta_K_quadrilatero, monta_F_quadrilatero

using BenchmarkTools: @benchmarkable, BenchmarkGroup
import BenchmarkTools as BT

const names = [
    :baseline, :ref, :iter, :iter_ref, :bacarmo,
]
const tests = [
    :build_vec, :build_mat, :build_both,
]
const small_tests_groups = (;
    build_small_vec  = [
        (:base        , build_small_vec_2d),
        (:base_precalc, build_small_vec_2d),
        (:ref         , build_small_vec_2d_ref!),
        (:ref_precalc , build_small_vec_2d_ref!),
        (:bacarmo     , monta_Fᵉ_quadrilatero!),
    ],
    build_small_mat  = [
        (:base        , build_small_mat_2d),
        (:base_precalc, build_small_mat_2d),
        (:ref         , build_small_mat_2d_ref!),
        (:ref_precalc , build_small_mat_2d_ref!),
        (:bacarmo     , monta_Kᵉ_quadrilatero!),
    ],
)
const tests_groups = (;
    build_vec  = [
        (:baseline, build_vec_2d),
        (:ref     , build_vec_2d_ref),
        (:iter    , build_vec_2d_iter),
        (:iter_ref, build_vec_2d_iterref),
    ],
    build_mat  = [
        (:baseline, build_mat_2d),
        (:ref     , build_mat_2d_ref),
        (:iter    , build_mat_2d_iter),
        (:iter_ref, build_mat_2d_iterref),
    ],
    build_both = [
        (:baseline, build_vec_mat_2d),
        (:ref     , build_vec_mat_2d_ref),
        (:iter    , build_vec_mat_2d_iter),
        (:iter_ref, build_vec_mat_2d_iterref),
    ],
)
const tests_tags = [
    (:build_small_vec,  ["build_small_vec"]),
    (:build_small_mat,  ["build_small_mat"]),
    (:build_vec,        ["build_vec"]),
    (:build_mat,        ["build_mat"]),
    (:build_both,       ["build_both"]),
]

const names_fullfuncs = (;
    baseline = build_vec_mat_2d,
    ref      = build_vec_mat_2d_ref,
    iter     = build_vec_mat_2d_iter,
    iter_ref = build_vec_mat_2d_iterref,
)
const names_tags = (;
    baseline = ["no_ref", "no_iter"],
    ref      = [   "ref", "no_iter"],
    iter     = ["no_ref",    "iter"],
    iter_ref = [   "ref",    "iter"],
    bacarmo  = ["bacarmo"          ],
)

suite = BenchmarkGroup()
for (test, test_tags) in tests_tags
    suite[test] = BenchmarkGroup(test_tags)
    if test in keys(tests_groups)
        for (impl_name, func) in tests_groups[test]
            suite[test][impl_name] = BenchmarkGroup(names_tags[impl_name])
        end
    end
end

const f = (x...) -> foldl(+, x)
const alpha, beta = 1.0, 1.0

const min_max = 2:5

begin
    local Ni = [1, 1]
    local hi = 1.0 ./ Ni

    local Xe, Ye = ([
        0.0, 1.0, 1.0, 0.0
    ], [
        0.0, 0.0, 1.0, 1.0
    ])

    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)

    local x2xis = x2xis_f(phis, Xe, Ye)
    local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

    for (name, func) in small_tests_groups[:build_small_vec]
        local sub_suite = suite[:build_small_vec]
        if name == :base
            sub_suite[name] = @benchmarkable begin
                local x2xis = x2xis_f($phis, $Xe, $Ye)
                local dx2xis = dx2xis_f($phi_derivs, $Xe, $Ye)

                $(func)(
                    $f,
                    x2xis, dx2xis,
                    $phis, $ws, $gauss_n
                )
            end
        elseif name == :base_precalc
            sub_suite[name] = @benchmarkable begin
                $(func)(
                    $f,
                    $x2xis, $dx2xis,
                    $phis, $ws, $gauss_n
                )
            end
        elseif name == :ref
            local ref_out = Ref(fill(0.0, (4,)))
            sub_suite[name] = @benchmarkable begin
                local x2xis = x2xis_f($phis, $Xe, $Ye)
                local dx2xis = dx2xis_f($phi_derivs, $Xe, $Ye)

                $(func)(
                    $ref_out,
                    $f,
                    x2xis, dx2xis,
                    $phis, $ws, $gauss_n
                )
                $ref_out[]
            end
        elseif name == :ref_precalc
            local ref_out = Ref(fill(0.0, (4,)))
            sub_suite[name] = @benchmarkable begin
                $(func)(
                    $ref_out,
                    $f,
                    $x2xis, $dx2xis,
                    $phis, $ws, $gauss_n
                )
                $ref_out[]
            end
        elseif name == :bacarmo
            local mat = fill(0.0, (4,))
            sub_suite[name] = @benchmarkable begin
                $(func)(
                    $mat,
                    $f,
                    $Xe, $Ye,
                    $ps, $ws
                )
                $mat
            end
        else
            error("Unknown symbol $impl")
        end
    end

    for (name, func) in small_tests_groups[:build_small_mat]
        local sub_suite = suite[:build_small_mat]
        if name == :base
            sub_suite[name] = @benchmarkable begin
                local dx2xis = dx2xis_f($phi_derivs, $Xe, $Ye)

                $(func)(
                    $alpha, $beta,
                    dx2xis,
                    $phis, $phi_derivs,
                    $ws, $gauss_n,
                )
            end
        elseif name == :base_precalc
            sub_suite[name] = @benchmarkable begin
                $(func)(
                    $alpha, $beta,
                    $dx2xis,
                    $phis, $phi_derivs,
                    $ws, $gauss_n,
                )
            end
        elseif name == :ref
            local ref_out = Ref(fill(0.0, (4,4)))
            sub_suite[name] = @benchmarkable begin
                local dx2xis = dx2xis_f($phi_derivs, $Xe, $Ye)

                $(func)(
                    $ref_out,
                    $alpha, $beta,
                    dx2xis,
                    $phis, $phi_derivs,
                    $ws, $gauss_n,
                )
                $ref_out
            end
        elseif name == :ref_precalc
            local ref_out = Ref(fill(0.0, (4,4)))
            sub_suite[name] = @benchmarkable begin
                $(func)(
                    $ref_out,
                    $alpha, $beta,
                    $dx2xis,
                    $phis, $phi_derivs,
                    $ws, $gauss_n,
                )
                $ref_out
            end
        elseif name == :bacarmo
            local mat = fill(0.0, (4,4))
            sub_suite[name] = @benchmarkable begin
                $(func)(
                    $mat,
                    $alpha, $beta,
                    $Xe, $Ye,
                    $ps, $ws,
                )
                $mat
            end
        else
            error("Unknown symbol $impl")
        end
    end
end

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
