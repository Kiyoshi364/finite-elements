include("finite-elements.jl")
include("../common.jl")

import .Common: calc_error, n_points_from_to

using .FiniteElements: fe_setup, fe_step
using .FiniteElements: example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = example(0x1)

min_max = 2:10

N_es = (1 .<< min_max)
hs = 1 ./ N_es

taus = hs
tss = (tau -> 0:tau:ex.T).(taus)
N_ts = Int.(floor.(ex.T ./ taus))

ABc01Ems = fe_setup.(ex, taus, hs, N_es)

errss = ((ABc01Em, args...) ->
((A, B, c_init0, c_init1, EQoLG, m, h, N_e, tau, ts, N_t) -> begin
    c_1 = c_init0
    c0 = c_init1
    errs = fill(0.0, (N_t+1,))

    errs[1] = Common.gauss_error(x -> exact(x, ts[1]), c_1, h)

    errs[2] = Common.gauss_error(x -> exact(x, ts[2]), c0, h)

    for i in 2:N_t
        t0 = ts[i]
        t1 = ts[i+1]
        c = fe_step(ex, A, B, c_1, c0, t0, tau, h, N_e, EQoLG, m)
        c_1, c0 = c0, c

        errs[i+1] = Common.gauss_error(x -> exact(x, t1), c, h)
    end
    maximum(errs)
end)(ABc01Em..., args...)).(ABc01Ems, hs, N_es, taus, tss, N_ts)

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
