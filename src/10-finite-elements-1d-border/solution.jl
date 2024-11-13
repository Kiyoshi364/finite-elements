include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: finite_elements, example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = example(0x1)

i = 3
N_e = (1 << i)
h = 1.0 / N_e
c = finite_elements(ex, h, N_e)

xs = Common.n_points_from_to(N_e-1, from=ex.x_begin, to=ex.x_end)
u = broadcast.(exact, xs)

err = Common.gauss_error(exact, c, h,
    ux_begin=ex.ux_begin, ux_end=ex.ux_end,
)

display(DataFrame(x=xs, solution=c, exact=u, diff=(u-c)))
println("\nError: ", err)

plotN = 1 << 8
many_xs = Common.n_points_from_to(plotN,
    from=ex.x_begin, to=ex.x_end,
    i_start=0, i_end=plotN
)
p = plot(many_xs, exact.(many_xs),
    label=L"$exact(x)$",
    xlabel=L"x",
    legend=:topleft,
)
plot!(p, cat(ex.x_begin, xs, ex.x_end, dims=1),
    cat(ex.ux_begin, c, ex.ux_end, dims=1),
    label=L"$solution(x)$",
    markershape=:circle,
)
savefig(p, "out.pdf")
