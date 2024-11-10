include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: finite_elements_setup
using .FiniteElements: generate_space, gauss_error_2d
using .FiniteElements: example

using Plots
using LaTeXStrings
using DataFrames

const exact, ex = example(0x0)

const noise = true
# const min_max = 2:7
const min_max = 2:6
const Nis = [ [ 1 << i for _ in 1:2 ] for i in min_max ]

const his = broadcast(./, 1.0, Nis)
const hs = (hi -> sqrt(foldl(+, hi .* hi))).(his)

const errss = @time (( (hi, Ni) -> begin
    local X, Y = generate_space(hi, Ni, noise=noise)
    local K, F, LG, EQoLG, m = finite_elements_setup(ex, X, Y, Ni)
    local c = K \ F
    gauss_error_2d(exact, c, X, Y, Ni, LG, EQoLG)
end).(his, Nis))

display(DataFrame(N=min_max, hi=his, error=errss))

const p = plot(
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
    label="finite element errors",
    markershape=:circle,
)
savefig(p, "out.pdf")
