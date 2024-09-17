include("finite-elements.jl")
include("../common.jl")

import .Common: calc_error, n_points_from_to

using .FiniteElements: fe_setup, fe_step
using .FiniteElements: bacarmo_example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = bacarmo_example(0x0)

min_max = 2:10

N_es = (1 .<< min_max)
N_ts = N_es
hs = 1 ./ N_es
taus = hs

ABcEs = fe_setup.(ex, taus, hs, N_es)

tss = (tau -> 0:tau:ex.T).(taus)

errss = ((ABcE, args...) ->
((A, B, c_init, EQoLG, h, tau, ts, N_e, N_t) -> begin
    c0 = c_init
    errs = fill(0.0, (N_t + 1,))

    errs[1] = Common.gauss_error(x -> exact(x, 0.0), c0, h)

    for i in 1:N_t
        t0 = ts[i]
        t1 = ts[i+1]
        c = fe_step(ex, A, B, c0, t0, tau, h, N_e, EQoLG)
        c0 = c

        errs[i+1] = Common.gauss_error(x -> exact(x, t1), c, h)
    end
    maximum(errs)
end)(ABcE..., args...)).(ABcEs, hs, taus, tss, N_es, N_ts)

display(DataFrame(N=min_max, h=hs, error=errss))

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
