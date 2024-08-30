include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: finite_elements, example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = example(0x0, 0x0)

i = 3
N = (1 << i) - 1
h = (ex.x_end - ex.x_begin) / (N + 1)
c = finite_elements(ex, h, N)

xs = Common.n_points_from_to(N, from=ex.x_begin, to=ex.x_end)
u = broadcast.(exact, xs)

err = Common.calc_error(c, u)

display(DataFrame(x=xs, solution=c, exact=u))

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
