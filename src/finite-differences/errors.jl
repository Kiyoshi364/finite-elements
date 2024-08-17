include("finite-difference.jl")

using .FiniteDifferences: finite_differences, example

using Plots
using LaTeXStrings
using DataFrames

function n_points_from_to(n; from=0, to=1)
    inter = (1:n) ./ (n+1)
    xs = (to - from) * inter .+ from
    xs
end

ex = example(0x0, 0x0)

min_max = 2:10

Ns = (1 .<< min_max) .- 1
hs = (ex.x_end - ex.x_begin) ./ (Ns .+ 1)
uhs = finite_differences.(ex, hs, Ns)

xss = n_points_from_to.(Ns, from=ex.x_begin, to=ex.x_end)
us = broadcast.(ex.exact, xss)

# errs = norm.(uhs .- us) ./ norm.(us)
errs = maximum.(broadcast.(abs, (uhs .- us)))

display(DataFrame(h=hs, error=errs))

p = plot(hs, hs .* hs,
    label=L"$O(h^2)$",
    yscale=:log10,
    xscale=:log10,
    xlabel=L"log_2(h)",
    legend=:topleft,
)
plot!(p, hs, errs,
    label="finite-difference errors",
    markershape=:circle,
)
savefig(p, "out.pdf")
