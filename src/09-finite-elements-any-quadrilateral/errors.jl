include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: build_vec_mat_2d, build_vec_mat_2d_ref
using .FiniteElements: build_vec_mat_2d_iter, build_vec_mat_2d_iterref
using .FiniteElements: generate_space, gauss_error_2d
using .FiniteElements: example

using .FiniteElements: gauss_quadrature_table
using .FiniteElements: phis_f, phi_derivs_f
using .FiniteElements: build_LG, build_EQ

include("bacarmo.jl")
using .Bacarmo: monta_F_quadrilatero, monta_K_quadrilatero

using Plots
using LaTeXStrings
using DataFrames

function wrap_build(to_wrap :: Function) :: Function
(
    f, alpha, beta,
    X, Y,
    N_e, LG, EQ, EQoLG, m,
    phis, phi_derivs,
    ws, gauss_n,
) -> begin
    to_wrap(
        f, alpha, beta,
        X, Y,
        N_e, LG,
        EQoLG, m,
        phis, phi_derivs,
        ws, gauss_n,
    )
end end

function build_vec_mat_2d_bacarmo(
    f, alpha, beta,
    X, Y,
    N_e, LG, EQ, EQoLG, m,
    phis, phi_derivs,
    ws, gauss_n,
)
    local K = monta_K_quadrilatero(
        alpha, beta,
        X, Y,
        m, EQ, LG,
    )
    local F = monta_F_quadrilatero(
        f,
        X, Y,
        m, EQ, LG,
    )
    F, K
end

const names_funcs = [
    ("baseline", wrap_build(build_vec_mat_2d)),
    ("ref"     , wrap_build(build_vec_mat_2d_ref)),
    ("iter"    , wrap_build(build_vec_mat_2d_iter)),
    ("iter_ref", wrap_build(build_vec_mat_2d_iterref)),
    # ("bacarmo",  build_vec_mat_2d_bacarmo),
]

const out_dir = "./out/"
const ext = "svg"

const exact, ex = example(0x0)

const noise = true
const min_max = 2:9
const Nis = [ [ 1 << i for _ in 1:2 ] for i in min_max ]

const his = broadcast(./, 1.0, Nis)
const hs = (hi -> sqrt(foldl(+, hi .* hi))).(his)

const gauss_n = 5

const ws, ps = gauss_quadrature_table[gauss_n]

const phis = phis_f(ps)
const phi_derivs = phi_derivs_f(ps)

for (name, func) in names_funcs
    local errss = @time (( (hi, Ni) -> begin

        local N_e = foldl(*, Ni)

        local LG = build_LG(Ni)

        local m, EQ = build_EQ(Ni)

        local EQoLG = EQ[LG]

        local X, Y = generate_space(hi, Ni, noise=noise)

        local F, K = func(
            ex.f, ex.alpha, ex.beta,
            X, Y,
            N_e, LG, EQ, EQoLG, m,
            phis, phi_derivs,
            ws, gauss_n,
        )

        local c = K \ F
        gauss_error_2d(exact, c, X, Y, Ni, LG, EQoLG)
    end).(his, Nis))

    display(DataFrame(N=min_max, hi=his, error=errss))

    local p = plot(
        yscale=:log2,
        xscale=:log2,
        xlabel=L"log_2(h)",
        ylabel=L"log_2",
        legend=:topleft,
    )
    plot!(p, hs, hs .* hs,
        label=L"$O(h^2)$",
    )
    plot!(p, hs, errss,
        label="max errors",
        markershape=:circle,
    )
    savefig(p, "$(out_dir)error-$(name).$(ext)")
end
