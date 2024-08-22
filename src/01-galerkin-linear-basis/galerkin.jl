module Galerkin

include("../examples.jl")
using .Examples: Example, example as common_example

include("../common.jl")
using .Common: build_tridiagonal, n_points_from_to

export galerkin, example

function build_mat(alpha, beta, h, dim)
    @assert dim >= 1
    @assert alpha > 0
    a = (- alpha / h) + (beta * h / 6)
    b = (2 * alpha / h) + (2 * beta * h / 3)
    build_tridiagonal(b, a, dim)
end

function build_vec(f, h, dim;
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert dim >= 2

    xs = n_points_from_to(dim, from=x_begin, to=x_end)
    f.(xs, h)
end

function galerkin(ex :: Example, h, N)
    galerkin(ex.f, ex.alpha, ex.beta, h, N;
        ux_begin=ex.ux_begin, ux_end=ex.ux_end,
        x_begin=ex.x_begin, x_end=ex.x_end,
    )
end

function galerkin(f, alpha, beta, h, N;
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert h == (x_end - x_begin) / (N + 1)
    K = build_mat(alpha, beta, h, N)
    F = build_vec(f, h, N,
        ux_begin=ux_begin, ux_end=ux_end,
        x_begin=x_begin, x_end=x_end,
    )
    c = K \ F
    c
end

function example(f_index :: UInt8, var_index :: UInt8 = 0) :: Tuple{Any, Example}
    func, ex = common_example(f_index, var_index)

    f = (
      f_index == 0 ? ((x, h) -> x * h) :
      f_index == 1 ? ((x, h) -> 8 * h) : # Note: only works for alpha=1 and beta=0
      error("function_index not supported")
    )

    func.exact, Example(ex,
        f=f
    )
end

end # module Galerkin
