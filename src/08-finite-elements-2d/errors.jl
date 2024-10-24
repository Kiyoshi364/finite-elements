include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: finite_elements_setup, gauss_error_2d
using .FiniteElements: example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = example(0x0)

min_max = 2:7
Nis = [ [ 1 << i for _ in 1:2 ] for i in min_max ]

N_es = foldl.(*, Nis)
his = broadcast(./, 1, Nis)
hs = (hi -> sqrt(foldl(+, hi .* hi))).(his)

KsFsEsms = finite_elements_setup.(ex, his, Nis)

cs = ( (KFEm) ->
((K, F, EQoLG, m) ->
    K \ F
)(KFEm...)).(KsFsEsms)

errss = ( (KFEm, args...) ->
((K, F, EQoLG, m, c, hi, Ni) ->
    gauss_error_2d(exact, c, hi, Ni, EQoLG)
)(KFEm..., args...)).(KsFsEsms, cs, his, Nis)

display(DataFrame(N=min_max, hi=his, error=errss))

p = plot(
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
