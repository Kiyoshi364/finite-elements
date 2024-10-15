include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: fe_setup, fe_step
using .FiniteElements: example

using Plots
using LaTeXStrings
using DataFrames

fake_exact, ex = example(0x2)

i = 4
N_e = (1 << i)
h = 1.0 / N_e
xs = Common.n_points_from_to(N_e-1)

tau = h
ts = 0:tau:ex.T
N_t = Int(floor(ex.T / tau))

A, B, c0, c1, EQoLG, m = fe_setup(ex, tau, h, N_e)

ylims=(0,0.11)

errs = fill(0.0, (N_t+1,))
c_1 = c0
anim = @animate for i in 0:N_t
    t1 = ts[i+1]

    c = ((i <= 0) ? begin
        c0
    end : (i == 1) ? begin
        c1
    end : begin
        t0 = ts[i]
        c = fe_step(ex, A, B, c_1, c0, t0, tau, h, N_e, EQoLG, m)

        c
    end)

    println("[$i/$N_t]")

    p = plot(
        legend=:topleft,
        xlabel=L"x",
        ylim=ylims,
    )
    plot!(p, cat(0.0, xs, 1.0, dims=1),
        cat(fake_exact(0.0, t1), c, fake_exact(1.0, t1), dims=1),
        label="solution(x, $t1)",
        markershape=:circle,
        color=2,
    )
    savefig(p, "out-$i.pdf")

    global c_1, c0 = c0, c
end

gif(anim, "out.gif", fps=1, show_msg=false)
