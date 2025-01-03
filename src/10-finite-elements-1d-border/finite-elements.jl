module FiniteElements

include("../examples.jl")
using .Examples: Example, example_begin_end as example, bacarmo_example

include("../common.jl")
using .Common: gauss_quadrature_table

using SparseArrays: spzeros

export finite_elements
export Example, example

phi = [ (xi -> (1 - xi) / 2), (xi -> (1 + xi) / 2) ]

phi_deriv = [ (xi -> - 0.5), (xi -> 0.5) ]

function build_small_mat(alpha, beta, h;
    gauss_n = 2,
)
    dim = 2

    ws, ps = gauss_quadrature_table[gauss_n]

    K_alpha = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p = ps[g_i]
                K_alpha[i,j] += ws[j] * (phi_deriv[i](p)*phi_deriv[j](p))
            end
        end
    end
    K_alpha *= 2 * alpha / h

    K_beta = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p = ps[g_i]
                K_beta[i,j] += ws[g_i] * (phi[i](p)*phi[j](p))
            end
        end
    end
    K_beta *= beta * h / 2

    K_alpha + K_beta
end

function build_mat(alpha, beta, h, N_e, EQoLG, m;
    gauss_n = 2,
)
    K_e = build_small_mat(alpha, beta, h,
        gauss_n=gauss_n,
    )

    K = spzeros((m+1, m+1))
    for e in 1:N_e
        _1 = EQoLG[1, e]
        _2 = EQoLG[2, e]
        K[_1,_1] += K_e[1,1]
        K[_1,_2] += K_e[1,2]
        K[_2,_1] += K_e[2,1]
        K[_2,_2] += K_e[2,2]
    end
    K[begin:end-1, begin:end-1]
end

function build_small_vec(f, h, e;
    gauss_n = 5,
)
    x2xi = (xi, e) -> h * ((1 + xi) / 2 + (e - 1))

    ws, ps = gauss_quadrature_table[gauss_n]
    F = fill(0.0, (2,))
    for i in 1:2
        for g_i in 1:gauss_n
            p = ps[g_i]
            F[i] += ws[g_i] * (f(x2xi(p, e))*phi[i](p))
        end
    end

    h/2 * F
end

function build_vec(f, h, N_e, EQoLG, m;
    gauss_n = 5,
)
    F = fill(0.0, (m+1,))
    for e in 1:N_e
        F_e = build_small_vec(f, h, e,
            gauss_n=gauss_n,
        )
        _1 = EQoLG[1, e]
        _2 = EQoLG[2, e]
        F[_1] += F_e[1]
        F[_2] += F_e[2]
    end

    F[begin:end-1]
end

function finite_elements(ex :: Example, h, N)
    finite_elements(ex.f, ex.alpha, ex.beta,
        ex.ux_begin, ex.ux_end, h, N)
end

function finite_elements(f, alpha, beta, p, q, h, N_e)

    LG = transpose(cat(1:N_e, 2:N_e+1, dims=2))
    EQ = cat(
        N_e,
        1:N_e-1,
        N_e,
        dims=1
    )
    m = N_e - 1
    EQoLG = EQ[LG]

    K = build_mat(alpha, beta, h, N_e, EQoLG, m)
    F = build_vec(f, h, N_e, EQoLG, m)

    K_e = build_small_mat(alpha, beta, h)
    F[begin] -= p * K_e[2, 1]
    F[end] -= q * K_e[1, 2]

    c = K \ F
    c
end

end # module FiniteElements
