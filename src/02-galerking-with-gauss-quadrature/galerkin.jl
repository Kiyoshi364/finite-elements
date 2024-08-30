module Galerkin

include("../examples.jl")
using .Examples: Example, example as common_example, bacarmo_example

include("../common.jl")
using .Common: build_tridiagonal, gauss_quadrature_table

export galerkin
export Example, example, bacarmo_example

function build_mat(alpha, beta, h, dim)
    @assert dim >= 1
    @assert alpha > 0
    a = (- alpha / h) + (beta * h / 6)
    b = (2 * alpha / h) + (2 * beta * h / 3)
    build_tridiagonal(b, a, dim)
end

function build_vec(f, h, dim;
    gauss_n = 5,
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert dim >= 2

    phi1 = xi -> (1 - xi) / 2
    phi2 = xi -> (1 + xi) / 2

    x2xi = (xi, i) -> (h * ((1 + xi) / 2 + (i - 1))) * (x_end - x_begin) + x_begin

    ws, ps = gauss_quadrature_table[gauss_n]
    F = fill(0.0, (dim,))
    for i in 1:dim
        for j in 1:gauss_n
            F[i] += ws[j] * (f(x2xi(ps[j], i))*phi2(ps[j]) + f(x2xi(ps[j], i+1))*phi1(ps[j]))
        end
    end

    h/2 * F
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

    func.exact, Example(ex,
        f = x -> ((- ex.alpha) * func.deriv_2(x)) + (ex.beta * func.exact(x))
    )
end

end # module Galerkin
