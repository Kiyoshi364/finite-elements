module FiniteDifferences

include("../examples.jl")

using .Examples: Example, example

export finite_differences
export Example, example

using LinearAlgebra

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
    ret
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
    ret
end

function finite_differences(ex :: Example, h, N)
    finite_differences(ex.f, ex.alpha, ex.beta, h, N;
        ux_begin=ex.ux_begin, ux_end=ex.ux_end,
        x_begin=ex.x_begin, x_end=ex.x_end,
    )
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
    uh
end

end # module FiniteDifferences
