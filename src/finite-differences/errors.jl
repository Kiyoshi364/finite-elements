include("finite-difference.jl")

using .FiniteDifferences: finite_differences

using Plots
using LaTeXStrings
using DataFrames

function n_points_from_to(n; from=0, to=1)
    inter = (1:n) ./ (n+1)
    xs = (to - from) * inter .+ from
    return xs
end

x_begin, x_end = 0, 1
alpha, beta = 1, 1

if true
    exact(x) = x * (x-1)
    deriv_2(x) = 2
else
    exact(x) = sin(pi * x)
    deriv_2(x) = - pi * pi * sin(pi * x)
end
func(x) = (deriv_2(x) * (- alpha)) + (beta * exact(x))
@assert abs(exact(0.0) - 0.0) < 1e-10
@assert abs(exact(1.0) - 0.0) < 1e-10

min_max = 2:10

Ns = (1 .<< min_max) .- 1
hs = (x_end - x_begin) ./ (Ns .+ 1)
uhs = finite_differences.(func, alpha, beta, hs, Ns,
    ux_begin=exact(x_begin),
    ux_end=exact(x_end),
    x_begin=x_begin,
    x_end=x_end,
)

us = broadcast.(
    exact,
    n_points_from_to.(Ns, from=x_begin, to=x_end)
)

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
