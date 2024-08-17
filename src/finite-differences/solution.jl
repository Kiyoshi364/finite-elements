include("finite-difference.jl")
include("../common.jl")

import .Common: n_points_from_to

using .FiniteDifferences: finite_differences, example

using Plots
using LaTeXStrings
using DataFrames

ex = example(0x0, 0x0)

i = 4
N = (1 << i) - 1
h = (ex.x_end - ex.x_begin) / (N + 1)
uh = finite_differences(ex, h, N)

xs = Common.n_points_from_to(N, from=ex.x_begin, to=ex.x_end)
u = broadcast.(ex.exact, xs)

err = Common.calc_error(uh, u)

display(DataFrame(x=xs, solution=uh, exact=u))

plotN = 1 << 8
many_xs = Common.n_points_from_to(plotN,
    from=ex.x_begin, to=ex.x_end,
    i_start=0, i_end=plotN
)
p = plot(many_xs, ex.exact.(many_xs),
    label=L"$exact(x)$",
    xlabel=L"x",
    legend=:topleft,
)
plot!(p, cat(ex.x_begin, xs, ex.x_end, dims=1),
    cat(ex.ux_begin, uh, ex.ux_end, dims=1),
    label=L"$solution(x)$",
    markershape=:circle,
    seriestype=:scatter,
)
savefig(p, "out.pdf")
