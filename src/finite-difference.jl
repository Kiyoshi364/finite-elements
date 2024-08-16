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

function build_vec(f, alpha, a, b, hsqr, dim)
    @assert dim >= 2

    offset = cat(
        [f(a) * alpha],
        fill(0, (dim-2,)),
        [f(b) * alpha],
    dims=1)

    inter = (1:dim) ./ (dim+1)
    xs = (b - a) * inter .+ a
    fxs = f.(xs)
    fs = fxs * hsqr
    sum = fs .+ offset

    ret = sum
    return ret
end

function finite_difference(f, alpha, beta, r_begin, r_end, h, N)
    @assert h == (r_end - r_begin) / (N + 1)
    hsqr = h * h
    A = build_mat(alpha, beta, hsqr, N)
    b = build_vec(f, alpha, r_begin, r_end, hsqr, N)
    uh = A \ b
    return uh
end

bgn, nd = 0, 1
alpha = 1
beta = 2

if true
    exact(x) = x * (x-1)
    deriv_2(x) = -2
else
    exact(x) = sin(pi * x)
    deriv_2(x) = pi * pi * sin(pi * x)
end
func(x) = (deriv_2(x) * alpha) + (beta * exact(x))

min_max = 3:10

Ns = (1 .<< min_max) .- 1
hs = (nd - bgn) ./ (Ns .+ 1)
uhs = finite_difference.(func, alpha, beta, bgn, nd, hs, Ns)

us = build_vec.(exact, 0, bgn, nd, 1, Ns)

# errs = norm.(uhs .- us) ./ norm.(us)
errs = maximum.(broadcast.(abs, (uhs .- us)))

display(DataFrame(h=hs, error=errs))

p = plot(hs, hs .* hs,
    label=L"$O(h^2)$",
    yscale=:log10,
    xscale=:log10,
    xlabel=L"log_2(h)",
    legend=:bottomright,
)
plot!(p, hs, errs,
    label="finite-difference errors",
    markershape=:circle,
)
savefig(p, "out.pdf")

end # module FiniteDifferences
