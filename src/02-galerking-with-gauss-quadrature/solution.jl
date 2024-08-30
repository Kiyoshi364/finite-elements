include("galerkin.jl")
include("../common.jl")

using .Galerkin: galerkin, bacarmo_example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = bacarmo_example(0x0)

i = 3
N = (1 << i) - 1
h = (ex.x_end - ex.x_begin) / (N + 1)
c = galerkin(ex, h, N)

xs = Common.n_points_from_to(N, from=ex.x_begin, to=ex.x_end)
u = broadcast.(exact, xs)

err = Common.gauss_error(exact, c, h,
    x_begin=ex.x_begin, x_end=ex.x_end,
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
