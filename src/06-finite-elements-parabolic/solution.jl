include("finite-elements.jl")
include("../common.jl")

using .FiniteElements: fe_setup, fe_step
using .FiniteElements: bacarmo_example

using Plots
using LaTeXStrings
using DataFrames

exact, ex = bacarmo_example(0x1)

i = 2
N_e = (1 << i)
N_t = N_e
h = 1.0 / N_e
tau = 1.0 / N_t
A, B, c0, EQoLG = fe_setup(ex, tau, h, N_e)

ts = 0:tau:ex.T

xs = Common.n_points_from_to(N_e-1)

plotN = 1 << 8
many_xs = Common.n_points_from_to(plotN,
    i_start=0, i_end=plotN
)

ylims=(0,0.11)

errs = fill(0.0, (N_t,))
anim = @animate for i in 0:N_t
    t1 = ts[i+1]

    c = ((i <= 0) ? begin
        c0
    end : begin
        t0 = ts[i]
        c = fe_step(ex, A, B, c0, t0, tau, h, N_e, EQoLG)
        u = broadcast.(x -> exact(x, t1), xs)

        println("\nTime($i): $t0")
        display(DataFrame(x=xs, solution=c, exact=u, diff=(u-c)))

        errs[i] = Common.gauss_error(x -> exact(x, t1), c, h)
        println("\nError: ", errs[i])
        c
    end)

    p = plot(
        legend=:topleft,
        xlabel=L"x",
        ylim=ylims,
    )
    plot!(p, many_xs, exact.(many_xs, t1),
        label="exact(x, $t1)",
    )
    plot!(p, cat(0.0, xs, 1.0, dims=1),
        cat(exact(0.0, t1), c, exact(1.0, t1), dims=1),
        label="solution(x, $t1)",
        markershape=:circle,
    )
    savefig(p, "out-$i.pdf")

    global c0 = c
end

gif(anim, "out.gif", fps=1, show_msg=false)
