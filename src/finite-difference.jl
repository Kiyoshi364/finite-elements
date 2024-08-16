module FiniteDifferences

# begin utils
broadcast(f, xs) = f.(xs)
# end utils

using LinearAlgebra
using Plots
using LaTeXStrings
using DataFrames

function build_mat(alpha, beta, hsqr, dim)
    #    (,~ $ ] }.@,@# 1 2 1 ,:@, 0 $~ 0>.-&2)
    @assert dim >= 2
    @assert alpha > 0
    a = - alpha
    b = (2 * alpha) + (beta * hsqr)
    zeros = fill(0, (dim-2,))
    prefix = [a, b, a]

    template = cat(prefix, zeros, dims=1)
    templates = fill(template, (dim,))
    flat = cat(templates..., dims=1)
    un_headed = flat[2:end-(dim-1)]
    ret = reshape(un_headed, (dim, dim))
    return ret
end

function build_vec(f, hsqr, dim;
    alpha=1,
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert dim >= 2

    offset = cat(
        [ux_begin * alpha],
        fill(0, (dim-2,)),
        [ux_end * alpha],
    dims=1)

    inter = (1:dim) ./ (dim+1)
    xs = (x_end - x_begin) * inter .+ x_begin
    fxs = f.(xs)
    fs = fxs * hsqr
    sum = fs + offset

    ret = sum
    return ret
end

function finite_differences(f, alpha, beta, h, N;
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert h == (x_end - x_begin) / (N + 1)
    hsqr = h * h
    A = build_mat(alpha, beta, hsqr, N)
    b = build_vec(f, hsqr, N,
        alpha=alpha,
        ux_begin=ux_begin, ux_end=ux_end,
        x_begin=x_begin, x_end=x_end,
    )
    uh = A \ b
    return uh
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

us = build_vec.(exact, 1, Ns,
    x_begin=x_begin,
    x_end=x_end,
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

end # module FiniteDifferences
