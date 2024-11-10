include("finite-elements.jl")

using .FiniteElements: finite_elements_setup
using .FiniteElements: generate_space, gauss_error_2d
using .FiniteElements: example

using Plots

const exact, ex = example(0x0)

const noise = true
const i = 2
const Ni = [ 1 << i for _ in 1:2 ]
const N_e = foldl(*, Ni)
const hi = 1.0 ./ Ni

const X, Y = generate_space(hi, Ni, noise=noise)

const K, F, LG, EQoLG, m = finite_elements_setup(ex, X, Y, Ni)

const c = K \ F

display(c)
display([
    exact([X[idx], Y[idx]])
    for i in 2:Ni[2]
    for j in 2:Ni[1]
    for idx in [(i-1) * (Ni[1]+1) + (j-1) + 1]
])

const err = gauss_error_2d(exact, c, X, Y, Ni, LG, EQoLG)
println("\nError: ", err)

# solution = instantiate_solution(c, phi, hi, Ni, EQoLG)

# p = plot(
#     seriestype=:surface,
# )
# plot!(p, 0:(hi[1]/4):1, 0:(hi[2]/4):1,
#     ((x,y) -> begin
#         display((x, y))
#         solution([x, y])
#     end),
#     label="solution",
# )

# plot!(p, many_xs, exact.(many_xs, t1),
#     label="exact(x, $t1)",
# )
# plot!(p, cat(0.0, xs, 1.0, dims=1),
#     cat(exact(0.0, t1), c, exact(1.0, t1), dims=1),
#     label="solution(x, $t1)",
#     markershape=:circle,
# )
# savefig(p, "out-$i.pdf")
