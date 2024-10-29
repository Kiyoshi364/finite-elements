module FiniteElements

include("../examples2d.jl")
using .Examples2d: Example, example

include("../common.jl")
using .Common: gauss_quadrature_table

using SparseArrays: spzeros

export finite_elements_setup
export gauss_error_2d
export Example, example

const phi = [
    ((xi1, xi2) -> (1 - xi1)*(1 - xi2) / 4),
    ((xi1, xi2) -> (1 + xi1)*(1 - xi2) / 4),
    ((xi1, xi2) -> (1 + xi1)*(1 + xi2) / 4),
    ((xi1, xi2) -> (1 - xi1)*(1 + xi2) / 4),
]

const phi_deriv = [
    [
        ((xi1, xi2) -> (- 0.25) * (1 - xi2)),
        ((xi1, xi2) -> (  0.25) * (1 - xi2)),
        ((xi1, xi2) -> (  0.25) * (1 + xi2)),
        ((xi1, xi2) -> (- 0.25) * (1 + xi2)),
    ], [
        ((xi1, xi2) -> (- 0.25) * (1 - xi1)),
        ((xi1, xi2) -> (- 0.25) * (1 + xi1)),
        ((xi1, xi2) -> (+ 0.25) * (1 + xi1)),
        ((xi1, xi2) -> (+ 0.25) * (1 - xi1)),
    ]
]

const hipe_x2xi = (hi, pe) -> (xi1, xi2) -> [
    (hi[1] * ((1 + xi1) / 2) + pe[1]),
    (hi[2] * ((1 + xi2) / 2) + pe[2]),
]

function build_small_mat_2d(alpha, beta, hi;
    gauss_n = 2,
)
    dim = 4

    ws, ps = gauss_quadrature_table[gauss_n]

    K_alpha1 = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    K_alpha1[i,j] += ws[g_i] * ws[g_j] *
                        (phi_deriv[1][j](p1, p2)*phi_deriv[1][i](p1, p2))
                end
            end
        end
    end
    K_alpha1 *= alpha * hi[2] / hi[1]

    K_alpha2 = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    K_alpha2[i,j] += ws[g_i] * ws[g_j] *
                        (phi_deriv[2][j](p1, p2)*phi_deriv[2][i](p1, p2))
                end
            end
        end
    end
    K_alpha2 *= alpha * hi[1] / hi[2]

    K_beta = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    K_beta[i,j] += ws[g_i] * ws[g_j] *
                        (phi[j](p1, p2)*phi[i](p1, p2))
                end
            end
        end
    end
    K_beta *= beta * hi[1] * hi[2] / 4

    K_alpha1 + K_alpha2 + K_beta
end

function build_mat_2d(alpha, beta, hi, N_e, EQoLG, m;
    gauss_n = 2
)
    K_e = build_small_mat_2d(alpha, beta, hi,
        gauss_n=gauss_n,
    )

    K = spzeros((m+1, m+1))
    for e in 1:N_e
        _1 = EQoLG[1, e]
        _2 = EQoLG[2, e]
        _3 = EQoLG[3, e]
        _4 = EQoLG[4, e]
        K[_1,_1] += K_e[1,1]
        K[_1,_2] += K_e[1,2]
        K[_1,_3] += K_e[1,3]
        K[_1,_4] += K_e[1,4]
        K[_2,_1] += K_e[2,1]
        K[_2,_2] += K_e[2,2]
        K[_2,_3] += K_e[2,3]
        K[_2,_4] += K_e[2,4]
        K[_3,_1] += K_e[3,1]
        K[_3,_2] += K_e[3,2]
        K[_3,_3] += K_e[3,3]
        K[_3,_4] += K_e[3,4]
        K[_4,_1] += K_e[4,1]
        K[_4,_2] += K_e[4,2]
        K[_4,_3] += K_e[4,3]
        K[_4,_4] += K_e[4,4]
    end
    K[begin:end-1, begin:end-1]
end

function build_small_vec_2d(f, hi, pe;
    gauss_n = 5,
)
    dim = 4

    x2xi = hipe_x2xi(hi, pe)

    ws, ps = gauss_quadrature_table[gauss_n]
    F = fill(0.0, (dim,))
    for i in 1:dim
        for g_i in 1:gauss_n
            p1 = ps[g_i]
            for g_j in 1:gauss_n
                p2 = ps[g_j]
                F[i] += ws[g_i] * ws[g_j] * (
                    f(x2xi(p1, p2))*phi[i](p1, p2)
                )
            end
        end
    end
    F *= hi[1] * hi[2] / 4

    F
end

function build_vec_2d(f, hi, Ni, N_e, EQoLG, m;
    gauss_n = 5,
)
    F = fill(0.0, (m+1,))
    for e in 1:N_e
        pe = [
            (mod(e-1, Ni[1]) * hi[1]),
            (div(e-1, Ni[1]) * hi[2]),
        ]

        F_e = build_small_vec_2d(f, hi, pe,
            gauss_n=gauss_n,
        )
        _1 = EQoLG[1, e]
        _2 = EQoLG[2, e]
        _3 = EQoLG[3, e]
        _4 = EQoLG[4, e]
        F[_1] += F_e[1]
        F[_2] += F_e[2]
        F[_3] += F_e[3]
        F[_4] += F_e[4]
    end

    F[begin:end-1]
end

function build_LG(Ni)
    Nx, Ny = Ni
    i = x -> x
    front4 = i.(0:(Ny-1)) .* (Nx + 1) .+ 1
    _0to3 = i.(0:(Nx-1))

    topline = cat(
        broadcast(.+, front4, (_0to3,))...,
    dims=1)

    LG = transpose([ 0 ;; 1 ;; Nx+2 ;; Nx+1 ] .+ topline)

    LG
end

function build_EQ(Ni)
    Nx, Ny = Ni[1], Ni[2]
    i = x -> x
    m = (Nx-1) * (Ny-1)

    small_mat = broadcast(
        .+,
        i.(0:(Ny-2)) .* (Nx-1),
        (i.(1:(Nx-1)),)
    )
    small_mat_ext = cat.(small_mat, m+1, m+1, dims=1)

    EQ = cat(
        fill(m+1, (Nx+2)),
        small_mat_ext...,
        fill(m+1, Nx), dims=1
    )

    (m, EQ)
end

function finite_elements_setup(ex :: Example, hi, Ni)
    finite_elements_setup(ex.f, ex.alpha, ex.beta, hi, Ni)
end

function finite_elements_setup(f, alpha, beta, hi, Ni)
    N_e = foldl(*, Ni)

    LG = build_LG(Ni)

    m, EQ = build_EQ(Ni)

    EQoLG = EQ[LG]

    K = build_mat_2d(alpha, beta, hi, N_e, EQoLG, m)
    F = build_vec_2d(f, hi, Ni, N_e, EQoLG, m)

    K, F, EQoLG, m
end

function gauss_error_2d(exact, coefs, hi, Ni, EQoLG;
    gauss_n = 5,
)
    N_e = foldl(*, Ni)

    ws, ps = gauss_quadrature_table[gauss_n]

    coefs_ext = cat(coefs, 0.0, dims=1)

    acc = 0.0

    for e in 1:N_e
        pe = [
            (mod(e-1, Ni[1]) * hi[1]),
            (div(e-1, Ni[1]) * hi[2]),
        ]

        x2xi = hipe_x2xi(hi, pe)

        for g_i = 1:gauss_n
            p1 = ps[g_i]
            for g_j in 1:gauss_n
                p2 = ps[g_j]
                diff = exact(x2xi(p1, p2))
                for i in 1:4
                    diff -= coefs_ext[EQoLG[i, e]] * phi[i](p1, p2)
                end
                acc += ws[g_i] * ws[g_j] * diff * diff
            end
        end
    end
    acc *= foldl(*, hi) / 4.0

    sqrt(acc)
end

end # module FiniteElements
