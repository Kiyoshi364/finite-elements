module FiniteDifferences

include("../examples.jl")
using .Examples: Example, example as common_example, bacarmo_example

include("../common.jl")
using .Common: build_tridiagonal, n_points_from_to

export finite_differences
export Example, example, bacarmo_example

using LinearAlgebra

function build_mat(alpha, beta, hsqr, dim)
    @assert dim >= 1
    @assert alpha > 0
    a = - alpha
    b = (2 * alpha) + (beta * hsqr)
    build_tridiagonal(b, a, dim)
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

    xs = n_points_from_to(dim, from=x_begin, to=x_end)
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

function example(f_index :: UInt8, var_index :: UInt8 = 0) :: Tuple{Any, Example}
    func, ex = common_example(f_index, var_index)
    func.exact, Example(ex,
        f = x -> ((- ex.alpha) * func.deriv_2(x)) + (ex.beta * func.exact(x))
    )
end

end # module FiniteDifferences
