include("finite-elements.jl")

using .FiniteElements: gauss_quadrature_table
using .FiniteElements: generate_space

using .FiniteElements: build_LG, build_EQ
using .FiniteElements: phis_f, phi_derivs_f, x2xis_f, dx2xis_f

using .FiniteElements: build_small_vec_2d, build_small_vec_2d_ref!
using .FiniteElements: build_vec_2d, build_vec_2d_ref
using .FiniteElements: build_vec_2d_iter, build_vec_2d_iterref

using .FiniteElements: build_small_mat_2d, build_small_mat_2d_ref!
using .FiniteElements: build_mat_2d, build_mat_2d_ref
using .FiniteElements: build_mat_2d_iter, build_mat_2d_iterref

using .FiniteElements: build_vec_mat_2d, build_vec_mat_2d_ref
using .FiniteElements: build_vec_mat_2d_iter, build_vec_mat_2d_iterref

using LinearAlgebra: dot, norm

using BenchmarkTools: @benchmarkable
import BenchmarkTools

function run_tests(name :: String, func :: Function, max_i :: UInt8;
    bench_strategy :: Symbol = :all
)
    println("# Running tests: $name")
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
    println()
end

function run_bench(b :: BenchmarkTools.Benchmark;
    bench :: Bool = true
)
    b.params.evals = (bench) ? b.params.evals : 1
    BenchmarkTools.tune!(b)
    BenchmarkTools.run_result(b)
end

function print_trial(trial :: BenchmarkTools.Trial;
    bench :: Bool = true
)
    if (false && bench)
        display(trial)
    else
        local trialmin = trial
        local trialallocs = BenchmarkTools.allocs(trialmin)
        println(
            "  ",
            BenchmarkTools.prettytime(BenchmarkTools.time(trialmin)),
            " (",
            trialallocs,
            " allocation",
            trialallocs == 1 ? "" : "s",
            ": ",
            BenchmarkTools.prettymemory(BenchmarkTools.memory(trialmin)),
            ")",
        )
    end
end

const default_tol = 1e-12
function test_check(i :: UInt8, expected, actual;
    tol=default_tol,
    msg=nothing,
)
    local err = norm((expected - actual)[begin:end])
    if msg == nothing
        println(". test $i: err is $err, tol is $tol")
    else
        println(". test $i ($msg): err is $err, tol is $tol")
    end
    if err >= tol
        println("==============================")
        println("expected:")
        display(expected)
        println("actual:")
        display(actual)
        println("==============================")
    end
end

const test_small_vec_2d_max = 0x02
function test_small_vec_2d_template(the_func :: Function, is_ref :: Bool) :: Function
(i :: UInt8; bench :: Bool = true) -> begin
    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)

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

    local phi_derivs = phi_derivs_f(ps)

    local x2xis = x2xis_f(phis, Xe, Ye)
    local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

    local f = (i == 0) ? (
        (x...) -> 4.0 / foldl(*, hi)
    ) : (i == 1) ? ( begin
        local prod_hi = foldl(*, hi)
        (x...) -> 16.0 * 9.0 * foldl(*, x) / prod_hi / prod_hi
    end ) : (i == 2) ? (
        (x...) -> foldl(+, x)
    ) : error("Index out of bounds")

    local b = (is_ref) ? begin
        local ref_out = Ref(fill(0.0, (4,)))
        @benchmarkable begin
            $(the_func)(
                $ref_out,
                $f, $x2xis,
                $dx2xis,
                $phis, $ws, $gauss_n
            )
            $ref_out[]
        end
    end : (@benchmarkable $(the_func)(
        $f, $x2xis,
        $dx2xis,
        $phis, $ws, $gauss_n
    ))
    local trial, ans = run_bench(b, bench=bench)
    test_check(i, expected, ans)
    print_trial(trial, bench=bench)

    ans
end end

const test_vec_2d_max = 0x02
function test_vec_2d_template(the_func :: Function) :: Function
(i :: UInt8; bench :: Bool = true) -> begin
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

    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)

    local X, Y = (i <= 1) ? (
        generate_space(hi, Ni)
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
        (x...) -> 4.0 / foldl(*, hi)
    ) : (i == 1) ? ( begin
        local prod_hi = foldl(*, hi)
        (x...) -> 16.0 * 9.0 * foldl(*, x) / prod_hi / prod_hi
    end ) : (i == 2) ? (
        (x...) -> foldl(+, x)
    ) : error("Index out of bounds")

    local b = @benchmarkable $(the_func)(
        $f, $X, $Y, $N_e, $LG, $EQoLG, $m,
        $phis, $phi_derivs,
        $ws, $gauss_n,
    )
    local trial, ans = run_bench(b, bench=bench)
    test_check(i, expected, ans,
        tol=(i == 2) ? 1e-7 : default_tol
    )
    print_trial(trial, bench=bench)

    ans
end end

const test_small_mat_2d_max = 0x02
function test_small_mat_2d_template(the_func :: Function, is_ref :: Bool) :: Function
(i :: UInt8; bench :: Bool = true) -> begin
    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local alpha, beta = (i == 0) ? (
        6.0, 0.0
    ) : (i == 1) ? (
        0.0, 576.0
    ) : (i == 2) ? (
        1.0, 1.0
    ) : error("Index out of bounds")

    local hi = [0.25, 0.25]
    local Xe, Ye = (i <= 1) ? (
        [ 0.0, hi[1], hi[1], 0.0 ],
        [ 0.0, 0.0, hi[2], hi[2] ],
    ) : (i == 2) ? (
        [ 0.0, 2.0, 3.0, 1.0 ],
        [ 0.0, 0.0, 1.0, 1.0 ],
    ) : error("Index out of bounds")

    local expected = (i == 0) ? ( [
        (4.0) (-1.0) (-2.0) (-1.0);
        (-1.0) (4.0) (-1.0) (-2.0);
        (-2.0) (-1.0) (4.0) (-1.0);
        (-1.0) (-2.0) (-1.0) (4.0);
    ] ) : (i == 1) ? ( [
        4.0 2.0 1.0 2.0;
        2.0 4.0 2.0 1.0;
        1.0 2.0 4.0 2.0;
        2.0 1.0 2.0 4.0;
    ] ) : (i == 2) ? ( [
        (0.7222222222222222) (0.1111111111111111) (0.0555555555555555) (-0.3888888888888888);
        (0.1111111111111111) (1.7222222222222222) (-0.3888888888888888) (-0.9444444444444444);
        (0.0555555555555555) (-0.3888888888888888) (0.7222222222222222) (0.1111111111111111);
        (-0.3888888888888888) (-0.9444444444444444) (0.1111111111111111) (1.7222222222222222);
    ] ) : error("Index out of bounds")

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)
    local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

    local b = (is_ref) ? begin
        local ref_out = Ref(fill(0.0, (4, 4)))
        @benchmarkable begin
            $(the_func)(
                $ref_out,
                $alpha, $beta,
                $dx2xis, $phis, $phi_derivs, $ws, $gauss_n
            )
            $ref_out[]
        end
    end : (@benchmarkable $(the_func)(
        $alpha, $beta,
        $dx2xis, $phis, $phi_derivs, $ws, $gauss_n
    ))
    local trial, ans = run_bench(b, bench=bench)
    test_check(i, expected, ans)
    print_trial(trial, bench=bench)

    ans
end end

const test_mat_2d_max = 0x01
function test_mat_2d_template(the_func :: Function) :: Function
(i :: UInt8; bench :: Bool = true) -> begin
    local gauss_n = 5

    local alpha, beta = 1.0, 1.0

    local Ni = [4, 3]
    local N_e = foldl(*, Ni)
    local hi = 1 ./ Ni

    local LG = build_LG(Ni)
    local m, EQ = build_EQ(Ni)
    local EQoLG = EQ[LG]

    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)

    local X, Y = (i == 0) ? (
        generate_space(hi, Ni)
    ) : (i == 1) ? ( [
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

    local expected = (i == 0) ? ( [
          2.81481    -0.62963     0.0        -0.0462963  -0.344907    0.0;
         -0.62963     2.81481    -0.62963    -0.344907   -0.0462963  -0.344907;
          0.0        -0.62963     2.81481     0.0        -0.344907   -0.0462963;
         -0.0462963  -0.344907    0.0         2.81481    -0.62963     0.0;
         -0.344907   -0.0462963  -0.344907   -0.62963     2.81481    -0.62963;
          0.0        -0.344907   -0.0462963   0.0        -0.62963     2.81481;
    ] ) : (i == 1) ? ( [
          2.78503   -0.710679    0.0        -0.203308  -0.246537    0.0;
         -0.710679   2.93963    -0.692196   -0.453569   0.0280001  -0.421416;
          0.0       -0.692196    2.86423     0.0       -0.321263    0.0453884;
         -0.203308  -0.453569    0.0         2.82204   -0.608672    0.0;
         -0.246537   0.0280001  -0.321263   -0.608672   2.87005    -0.885907;
          0.0       -0.421416    0.0453884   0.0       -0.885907    2.91071;
    ] ) : error("Index out of bounds")

    local b = @benchmarkable $(the_func)(
        $alpha, $beta, $X, $Y, $N_e, $LG, $EQoLG, $m,
        $phis, $phi_derivs,
        $ws, $gauss_n,
    )
    local trial, ans = run_bench(b, bench=bench)
    test_check(i, expected, ans,
        tol=(i == 1) ? 1e-2 : 1e-4
    )
    print_trial(trial, bench=bench)

    ans
end end

const test_vec_mat_2d_max = 0x00
function test_vec_mat_2d_template(the_func :: Function) :: Function
(i :: UInt8; bench :: Bool = true) -> begin
    local gauss_n = 5

    local alpha, beta = 1.0, 1.0

    local Ni = [4, 3]
    local N_e = foldl(*, Ni)
    local hi = 1 ./ Ni

    local LG = build_LG(Ni)
    local m, EQ = build_EQ(Ni)
    local EQoLG = EQ[LG]

    local gauss_n = 5
    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)

    local X, Y = ( [
        0.0, 0.25, 0.5, 0.75, 1.0,
        0.0, 0.266168, 0.493792, 0.747176, 1.0,
        0.0, 0.275391, 0.521668, 0.708237, 1.0,
        0.0, 0.25, 0.5, 0.75, 1.0,
    ], [
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.333333, 0.352246, 0.36139, 0.326172, 0.333333,
        0.666667, 0.633228, 0.693524, 0.689905, 0.666667,
        1.0, 1.0, 1.0, 1.0, 1.0,
    ] )

    local expected_vec = [
        0.047603316256427886,
        0.0685848331877617,
        0.09299899556951925,
        0.07729379019471166,
        0.08630556638050947,
        0.11514211312826136
    ]
    local expected_mat = [
          2.78503   -0.710679    0.0        -0.203308  -0.246537    0.0;
         -0.710679   2.93963    -0.692196   -0.453569   0.0280001  -0.421416;
          0.0       -0.692196    2.86423     0.0       -0.321263    0.0453884;
         -0.203308  -0.453569    0.0         2.82204   -0.608672    0.0;
         -0.246537   0.0280001  -0.321263   -0.608672   2.87005    -0.885907;
          0.0       -0.421416    0.0453884   0.0       -0.885907    2.91071;
    ]
    local expected = (expected_vec, expected_mat)

    local f = (x...) -> foldl(+, x)

    local b = @benchmarkable $(the_func)(
        $f, $alpha, $beta,
        $X, $Y, $N_e, $LG, $EQoLG, $m,
        $phis, $phi_derivs,
        $ws, $gauss_n,
    )
    local trial, ans = run_bench(b, bench=bench)
    test_check(i, expected[1], ans[1],
        tol=1e-7,
        msg="vec",
    )

    test_check(i, expected[2], ans[2],
        tol=1e-2,
        msg="mat",
    )
    print_trial(trial, bench=bench)

    ans
end end

const names_localvec_isref = [
    ("test_small_vec_2d", build_small_vec_2d, false),
    ("test_small_vec_2d_ref", build_small_vec_2d_ref!, true),
]
for (name, func, is_ref) in names_localvec_isref
    run_tests(name, test_small_vec_2d_template(func, is_ref), test_small_vec_2d_max, bench_strategy=:all)
end

const names_funcvec = [
    ("test_vec_2d", build_vec_2d),
    ("test_vec_2d_ref", build_vec_2d_ref),
    ("test_vec_2d_iter", build_vec_2d_iter),
    ("test_vec_2d_iterref", build_vec_2d_iterref),
]

for (name, func) in names_funcvec
    run_tests(name, test_vec_2d_template(func), test_vec_2d_max, bench_strategy=:last)
end

const names_localmat_isref = [
    ("test_small_mat_2d", build_small_mat_2d, false),
    ("test_small_mat_2d_ref", build_small_mat_2d_ref!, true),
]
for (name, func, is_ref) in names_localmat_isref
    run_tests(name, test_small_mat_2d_template(func, is_ref), test_small_mat_2d_max, bench_strategy=:all)
end

const names_funcmat = [
    ("test_mat_2d", build_mat_2d),
    ("test_mat_2d_ref", build_mat_2d_ref),
    ("test_mat_2d_iter", build_mat_2d_iter),
    ("test_mat_2d_iterref", build_mat_2d_iterref),
]
for (name, func) in names_funcmat
    run_tests(name, test_mat_2d_template(func), test_mat_2d_max, bench_strategy=:last)
end

const names_funcboth = [
    ("test_vec_mat_2d", build_vec_mat_2d),
    ("test_vec_mat_2d_ref", build_vec_mat_2d_ref),
    ("test_vec_mat_2d_iter", build_vec_mat_2d_iter),
    ("test_vec_mat_2d_iterref", build_vec_mat_2d_iterref),
]
for (name, func) in names_funcboth
    run_tests(name, test_vec_mat_2d_template(func), test_vec_mat_2d_max, bench_strategy=:all)
end
