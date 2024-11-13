include("finite-elements.jl")
include("../common.jl")

import .Common: calc_error, n_points_from_to

using .FiniteElements: finite_elements, example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = example(0x1)

min_max = 2:10

N_es = (1 .<< min_max)
hs = 1.0 ./ N_es
cs = finite_elements.(ex, hs, N_es)

errs = Common.gauss_error.(exact, cs, hs,
    ux_begin=ex.ux_begin, ux_end=ex.ux_end,
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
