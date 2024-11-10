include("finite-elements.jl")

using .FiniteElements: gauss_quadrature_table
using .FiniteElements: phi, phi_deriv

using .FiniteElements: build_LG, build_EQ
using .FiniteElements: build_small_vec_2d, build_vec_2d

using LinearAlgebra: dot, norm

using BenchmarkTools: @btime

function run_tests(name, func, max_i :: UInt8;
    bench_strategy=:all
)
    println("Running tests: $name")
    for i in 0x0:max_i
        local bench = (bench_strategy == :all) ? (
            true
        ) : (bench_strategy == :none) ? (
            false
        ) : (bench_strategy == :first) ? (
            i == 0
        ) : (bench_strategy == :last) ? (
            i == max_i
        ) : error("Invalid bench_strategy")
        func(i, bench=bench)
    end
end

const default_tol = 1e-12
function test_check(i :: UInt8, expected, actual;
    tol=default_tol
)
    local err = norm((expected - actual)[begin:end])
    if err >= tol
        println("==============================")
        println("test $i: err is $err, tol is $tol")
        println("expected:")
        display(expected)
        println("actual:")
        display(actual)
        println("==============================")
    else
        println("test $i: err is $err, tol is $tol")
    end
end

const test_small_vec_2d_max = 0x02
function test_small_vec_2d(i :: UInt8;
    bench=true
)
    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local app = (f, xs...) -> f(xs...)
    local sdim = 2
    local dim = 4

    local hi = [0.25, 0.25]
    local Xe, Ye = (i <= 1) ? (
        [ 0.0, hi[1], hi[1], 0.0 ],
        [ 0.0, 0.0, hi[2], hi[2] ],
    ) : (i == 2) ? (
        [ 0.0, 2.0, 3.0, 1.0 ],
        [ 0.0, 0.0, 1.0, 1.0 ],
    ) : error("Index out of bounds")

    local expected = (i == 0) ? (
        [ 1.0, 1.0, 1.0, 1.0 ]
    ) : (i == 1) ? (
        [ 4.0, 8.0, 16.0, 8.0 ]
    ) : (i == 2) ? (
        [ 0.6666666666666666, 1.0, 1.3333333333333333, 1.0 ]
    ) : error("Index out of bounds")

    local x2xis = (Xe, Ye) -> [
        [
            dot.([Xe, Ye], (app.(phi, ([pi, pj],)),))
            for pj in ps
        ]
        for pi in ps
    ]

    local dx2xis = (Xe, Ye) -> [
        [
            [
                dot([Xe, Ye][a], app.(phi_deriv[b], ([pi, pj],)))
                for a in 1:sdim
                for b in 1:sdim
            ]
            for pj in ps
        ]
        for pi in ps
    ]

    local f = (i == 0) ? (
        x -> 4.0 / foldl(*, hi)
    ) : (i == 1) ? (
        x -> 16. * 9.0 * foldl(*, x) / foldl(*, hi) / foldl(*, hi)
    ) : (i == 2) ? (
        x -> foldl(+, x)
    ) : error("Index out of bounds")

    local ans = (bench) ? (@btime build_small_vec_2d(
        $f, $(x2xis(Xe, Ye)),
        $(dx2xis(Xe, Ye)), $ws, $ps, $gauss_n
    )) : (build_small_vec_2d(
        f, x2xis(Xe, Ye),
        dx2xis(Xe, Ye), ws, ps, gauss_n
    ))

    test_check(i, expected, ans)

    ans
end

const test_vec_2d_max = 0x02
function test_vec_2d(i :: UInt8; bench=true)
    local gauss_n = 5

    local Ni = (i == 0) ? (
        [4, 3]
    ) : (i == 1) ? (
        [4, 4]
    ) : (i == 2) ? (
        [4, 3]
    ) : error("Index out of bounds")
    local N_e = foldl(*, Ni)
    local hi = 1 ./ Ni

    local LG = build_LG(Ni)
    local m, EQ = build_EQ(Ni)
    local EQoLG = EQ[LG]

    local X, Y = (i <= 1) ? (
        repeat(0.0:(hi[1]):1.0, Ni[2]+1),
        cat(
            ((x->repeat(x:x, Ni[1]+1)).(0.0:(hi[2]):1.0))...,
            dims=1
        ),
    ) : (i == 2) ? ( [
        0.0, 0.25, 0.5, 0.75, 1.0,
        0.0, 0.266168, 0.493792, 0.747176, 1.0,
        0.0, 0.275391, 0.521668, 0.708237, 1.0,
        0.0, 0.25, 0.5, 0.75, 1.0,
    ], [
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.333333, 0.352246, 0.36139, 0.326172, 0.333333,
        0.666667, 0.633228, 0.693524, 0.689905, 0.666667,
        1.0, 1.0, 1.0, 1.0, 1.0,
    ] ) : error("Index out of bounds")

    local expected = (i == 0) ? (
        [ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0 ]
    ) : (i == 1) ? (
        [ 144.0, 288.0, 432.0, 288.0, 576.0,
          864.0, 432.0, 864.0, 1296.0 ]
    ) : (i == 2) ? ( [
        0.047603316256427886,
        0.0685848331877617,
        0.09299899556951925,
        0.07729379019471166,
        0.08630556638050947,
        0.11514211312826136
    ] ) : error("Index out of bounds")

    local f = (i == 0) ? (
        x -> 4.0 / foldl(*, hi)
    ) : (i == 1) ? (
        x -> 16. * 9.0 * foldl(*, x) / foldl(*, hi) / foldl(*, hi)
    ) : (i == 2) ? (
        x -> foldl(+, x)
    ) : error("Index out of bounds")

    local ans = (bench) ? (@btime build_vec_2d(
        $f, $X, $Y, $N_e, $LG, $EQoLG, $m,
        gauss_n=$gauss_n
    )) : (build_vec_2d(
        f, X, Y, N_e, LG, EQoLG, m,
        gauss_n=gauss_n
    ))

    test_check(i, expected, ans,
        tol=(i == 2) ? 1e-7 : default_tol
    )

    ans
end

run_tests("test_small_vec_2d", test_small_vec_2d, test_small_vec_2d_max, bench_strategy=:all)
run_tests("test_vec_2d", test_vec_2d, test_vec_2d_max, bench_strategy=:last)
