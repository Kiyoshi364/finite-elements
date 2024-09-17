module FiniteElements

include("../examples.jl")
using .TimeExamples: Example, bacarmo_example

include("../common.jl")
using .Common: gauss_quadrature_table, n_points_from_to

using LinearAlgebra: lu
using SparseArrays: spzeros

export finite_elements
export fe_setup, fe_setup_AB, fe_setup_c0, fe_step
export Example, example, bacarmo_example

phi = [ (xi -> (1 - xi) / 2), (xi -> (1 + xi) / 2) ]

phi_deriv = [ (xi -> - 0.5), (xi -> 0.5) ]

function build_small_mat(alpha, beta, gamma, h;
    gauss_n = 2,
)
    dim = 2

    ws, ps = gauss_quadrature_table[gauss_n]

    K_alpha = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p = ps[g_i]
                K_alpha[i,j] += ws[g_i] * (phi_deriv[i](p)*phi_deriv[j](p))
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

    K_gamma = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p = ps[g_i]
                K_gamma[i,j] += ws[g_i] * (phi_deriv[j](p)*phi[i](p))
            end
        end
    end
    K_gamma *= gamma

    K_alpha + K_beta + K_gamma
end

function build_mat(alpha, beta, gamma, h, N_e, EQoLG;
    gauss_n = 2
)
    K_e = build_small_mat(alpha, beta, gamma, h,
        gauss_n=gauss_n,
    )

    K = spzeros((N_e, N_e))
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
    x2xi = (xi, e) -> (h * ((1 + xi) / 2 + (e - 1)))

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

function build_vec(f, h, N_e, EQoLG;
    gauss_n = 5,
)
    F = fill(0.0, (N_e,))
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

function finite_elements(ex :: Example, tau, h, N_e)
    finite_elements(ex.f, ex.u0, ex.T, tau, ex.alpha, ex.beta, ex.gamma, h, N_e)
end

function finite_elements(f, u0, T, tau, alpha, beta, gamma, h, N_e)
    A, B, c0, EQoLG = fe_setup(f, u0, tau, alpha, beta, gamma, h, N_e)

    ts = 0:tau:T

    for t0 in ts[begin:end-1]
        c = fe_step(f, A, B, c0, t0, tau, h, N_e, EQoLG)
        c0 = c
    end
    c0
end

function fe_setup(ex :: Example, tau, h, N_e)
    fe_setup(ex.f, ex.u0, ex.alpha, ex.beta, ex.gamma, tau, h, N_e)
end

function fe_setup(f, u0, alpha, beta, gamma, tau, h, N_e)
    LG = transpose(cat(1:N_e, 2:N_e+1, dims=2))
    EQ = cat(
        N_e,
        1:N_e-1,
        N_e,
        dims=1
    )
    EQoLG = EQ[LG]

    A, B = fe_setup_AB(alpha, beta, gamma, tau, h, N_e, EQoLG)

    c0 = fe_setup_c0(u0, h, N_e, EQ)

    A, B, c0, EQoLG
end

function fe_setup_AB(alpha, beta, gamma, tau, h, N_e, EQoLG)

    M = build_mat(0.0, 1.0, 0.0, h, N_e, EQoLG)
    K = build_mat(alpha, beta, gamma, h, N_e, EQoLG)

    K_half_tau = (tau / 2) * K

    A0 = M + K_half_tau
    B0 = M - K_half_tau

    A = lu(A0)
    B = B0

    A, B
end

function fe_setup_c0(u0, h, N_e, EQ)
    xs = n_points_from_to(N_e-1,
        i_start=0, i_end=N_e
    )

    eq_xs = fill(0.0, (N_e,))
    for i in 1:(N_e+1)
        _i = EQ[i]
        eq_xs[_i] += xs[i]
    end
    eq_xs = eq_xs[begin:end-1]

    u0eq_xs = u0.(eq_xs)
    u0eq_xs
end

function fe_step(ex :: Example, A, B, c0, t0, tau, h, N_e, EQoLG)
    fe_step(ex.f, A, B, c0, t0, tau, h, N_e, EQoLG)
end

function fe_step(f, A, B, c0, t0, tau, h, N_e, EQoLG)
    @assert h == 1 / N_e

    t0_half = t0 + (tau / 2)

    F = build_vec(x -> f(x, t0_half), h, N_e, EQoLG)

    c = A \ ((B * c0) + (tau * F))
    c
end

end # module FiniteElements
