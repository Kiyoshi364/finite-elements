include("finite-elements.jl")

using .FiniteElements: finite_elements_setup, gauss_error_2d
using .FiniteElements: example

using Plots

exact, ex = example(0x0)

i = 2
Ni = [ 1 << i for _ in 1:2 ]
N_e = foldl(*, Ni)
hi = 1.0 ./ Ni

K, F, EQoLG, m = finite_elements_setup(ex, hi, Ni)

c = K \ F

display(c)
display([
    exact([x1, x2])
    for x2 in hi[2]:hi[2]:1-hi[2]
    , x1 in hi[1]:hi[1]:1-hi[1]
])

err = gauss_error_2d(exact, c, hi, Ni, EQoLG)
println("\nError: ", err)

solution = instantiate_solution(c, phi, hi, Ni, EQoLG)
