include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: finite_elements, bacarmo_example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = bacarmo_example(0x0)

i = 3
N_e = (1 << i)
h = 1.0 / N_e
tau = h
c = finite_elements(ex, tau, h, N_e)

# TODO
xs = Common.n_points_from_to(N_e-1)
u = broadcast.(x -> exact(x, ex.T), xs)

err = Common.gauss_error(x -> exact(x, ex.T), c, h)

display(DataFrame(x=xs, solution=c, exact=u, diff=(u-c)))
println("\nError: ", err)

plotN = 1 << 8
many_xs = Common.n_points_from_to(plotN,
    i_start=0, i_end=plotN
)
p = plot(many_xs, exact.(many_xs, ex.T),
    label=L"$exact(x)$",
    xlabel=L"x",
    legend=:topleft,
)
plot!(p, cat(0.0, xs, 1.0, dims=1),
    cat(exact(0.0, ex.T), c, exact(1.0, ex.T), dims=1),
    label=L"$solution(x)$",
    markershape=:circle,
)
savefig(p, "out.pdf")
