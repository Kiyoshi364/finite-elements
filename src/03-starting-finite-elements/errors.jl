include("finite-elements.jl")
include("../common.jl")

import .Common: calc_error, n_points_from_to

using .FiniteElements: finite_elements, example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = example(0x0, 0x0)

min_max = 2:10

Ns = (1 .<< min_max) .- 1
hs = (ex.x_end - ex.x_begin) ./ (Ns .+ 1)
cs = finite_elements.(ex, hs, Ns)

errs = Common.gauss_error.(exact, cs, hs,
    x_begin=ex.x_begin,
    x_end=ex.x_end,
)

display(DataFrame(N=min_max, h=hs, error=errs))

p = plot(hs, hs .* hs,
    label=L"$O(h^2)$",
    yscale=:log2,
    xscale=:log2,
    xlabel=L"log_2(h)",
    ylabel=L"log_2",
    legend=:topleft,
)
plot!(p, hs, errs,
    label="finite element errors",
    markershape=:circle,
)
savefig(p, "out.pdf")
