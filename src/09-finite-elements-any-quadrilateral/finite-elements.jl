module FiniteElements

include("../examples2d.jl")
using .Examples2d: Example, example

include("../common.jl")
using .Common: gauss_quadrature_table

using LinearAlgebra: dot

using SparseArrays: spzeros

export finite_elements_setup
export generate_space
export gauss_error_2d
export phi, phi_deriv, instantiate_solution
export Example, example

# For tests
export build_LG, build_EQ
export x2xis_f, dx2xis_f
export build_small_vec_2d, build_vec_2d
export build_small_mat_2d, build_mat_2d

const phi = [
    (xi -> (1 - xi[1])*(1 - xi[2]) / 4),
    (xi -> (1 + xi[1])*(1 - xi[2]) / 4),
    (xi -> (1 + xi[1])*(1 + xi[2]) / 4),
    (xi -> (1 - xi[1])*(1 + xi[2]) / 4),
]

const phi_deriv = [
    [
        (xi -> (- 0.25) * (1 - xi[2])),
        (xi -> (  0.25) * (1 - xi[2])),
        (xi -> (  0.25) * (1 + xi[2])),
        (xi -> (- 0.25) * (1 + xi[2])),
    ], [
        (xi -> (- 0.25) * (1 - xi[1])),
        (xi -> (- 0.25) * (1 + xi[1])),
        (xi -> (+ 0.25) * (1 + xi[1])),
        (xi -> (+ 0.25) * (1 - xi[1])),
    ]
]

const sdim = 2
const app = (f, xs...) -> f(xs...)

const x2xis_f = (ps, Xe, Ye) -> [
    [
        dot.([Xe, Ye], (app.(phi, ([pi, pj],)),))
        for pj in ps
    ]
    for pi in ps
]

const dx2xis_f = (ps, Xe, Ye) -> [
    [
        [
            dot([Xe, Ye][a], app.(phi_deriv[b], ([pi, pj],)))
            for a in 1:sdim
            for b in 1:sdim
        ]
        for pj in ps
    ]
    for pi in ps
]

function build_small_mat_2d(alpha, beta,
    dx2xis, ws, ps, gauss_n,
)
    dim = 4

    K_alpha1 = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    _J = dx2xis[g_i][g_j]
                    J = 1.0 / ((_J[1] * _J[4]) - (_J[2] * _J[3]))
                    H_1 = (_J[4]*_J[4] + _J[2]*_J[2])
                    H_2 = - (_J[4]*_J[3] + _J[2]*_J[1])
                    K_alpha1[i,j] += J * ws[g_i] * ws[g_j] * (
                        phi_deriv[1][j]([p1, p2]) * (
                            (H_1 * phi_deriv[1][i]([p1, p2]))
                            + (H_2 * phi_deriv[2][i]([p1, p2]))
                        )
                    )
                end
            end
        end
    end
    K_alpha1 *= alpha

    K_alpha2 = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    _J = dx2xis[g_i][g_j]
                    J = 1.0 / ((_J[1] * _J[4]) - (_J[2] * _J[3]))
                    H_1 = - (_J[4]*_J[3] + _J[2]*_J[1])
                    H_2 = (_J[3]*_J[3] + _J[1]*_J[1])
                    K_alpha2[i,j] += J * ws[g_i] * ws[g_j] * (
                        phi_deriv[2][j]([p1, p2]) * (
                            (H_1 * phi_deriv[1][i]([p1, p2]))
                            + (H_2 * phi_deriv[2][i]([p1, p2]))
                        )
                    )
                end
            end
        end
    end
    K_alpha2 *= alpha

    K_beta = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    _J = dx2xis[g_i][g_j]
                    J = (_J[1] * _J[4]) - (_J[2] * _J[3])
                    K_beta[i,j] += J * ws[g_i] * ws[g_j] *
                        (phi[j]([p1, p2])*phi[i]([p1, p2]))
                end
            end
        end
    end
    K_beta *= beta

    K_alpha1 + K_alpha2 + K_beta
end

function build_mat_2d(alpha, beta, X, Y, N_e, LG, EQoLG, m;
    gauss_n = 5
)
    local sdim = 2
    local ws, ps = gauss_quadrature_table[gauss_n]

    local K = spzeros((m+1, m+1))
    for e in 1:N_e
        local LGe = LG[:, e]
        local Xe = X[LGe]
        local Ye = Y[LGe]

        local dx2xis = dx2xis_f(ps, Xe, Ye)

        local K_e = build_small_mat_2d(alpha, beta,
            dx2xis, ws, ps, gauss_n,
        )

        local _1 = EQoLG[1, e]
        local _2 = EQoLG[2, e]
        local _3 = EQoLG[3, e]
        local _4 = EQoLG[4, e]
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

function build_small_vec_2d(f, x2xis,
    dx2xis, ws, ps, gauss_n
)
    local dim = 4

    local F = fill(0.0, (dim,))
    for i in 1:dim
        for g_i in 1:gauss_n
            local p1 = ps[g_i]
            for g_j in 1:gauss_n
                local p2 = ps[g_j]
                local _x = x2xis[g_i][g_j]
                local _J = dx2xis[g_i][g_j]
                local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
                F[i] += J * ws[g_i] * ws[g_j] * (
                    f(_x)*phi[i]([p1, p2])
                )
            end
        end
    end

    F
end

function build_vec_2d(f, X, Y, N_e, LG, EQoLG, m;
    gauss_n = 5,
)
    local ws, ps = gauss_quadrature_table[gauss_n]

    local F = fill(0.0, (m+1,))
    for e in 1:N_e
        local LGe = LG[:, e]
        local Xe = X[LGe]
        local Ye = Y[LGe]

        local x2xis = x2xis_f(ps, Xe, Ye)
        local dx2xis = dx2xis_f(ps, Xe, Ye)

        local F_e = build_small_vec_2d(f, x2xis,
            dx2xis, ws, ps, gauss_n
        )
        local _1 = EQoLG[1, e]
        local _2 = EQoLG[2, e]
        local _3 = EQoLG[3, e]
        local _4 = EQoLG[4, e]
        F[_1] += F_e[1]
        F[_2] += F_e[2]
        F[_3] += F_e[3]
        F[_4] += F_e[4]
    end

    F[begin:end-1]
end

function build_LG(Ni)
    local Nx, Ny = Ni
    local i = x -> x
    local front4 = i.(0:(Ny-1)) .* (Nx + 1) .+ 1
    local _0to3 = i.(0:(Nx-1))

    local topline = cat(
        broadcast(.+, front4, (_0to3,))...,
    dims=1)

    local LG = transpose([ 0 ;; 1 ;; Nx+2 ;; Nx+1 ] .+ topline)

    LG
end

function build_EQ(Ni)
    local Nx, Ny = Ni[1], Ni[2]
    local i = x -> x
    local m = (Nx-1) * (Ny-1)

    local small_mat = broadcast(
        .+,
        i.(0:(Ny-2)) .* (Nx-1),
        (i.(1:(Nx-1)),)
    )
    local small_mat_ext = cat.(small_mat, m+1, m+1, dims=1)

    local EQ = cat(
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
    local N_e = foldl(*, Ni)

    local LG = build_LG(Ni)

    local m, EQ = build_EQ(Ni)

    local EQoLG = EQ[LG]

    local K = build_mat_2d(alpha, beta, hi, N_e, EQoLG, m)
    local F = build_vec_2d(f, hi, Ni, N_e, EQoLG, m)

    K, F, EQoLG, m
end

function generate_space(hi, Ni;
    noise=false
)
    X = repeat(0.0:(hi[1]):1.0, Ni[2]+1)
    Y = cat(
        ((x->repeat(x:x, Ni[1]+1)).(0.0:(hi[2]):1.0))...,
        dims=1
    )
    X, Y
end

function gauss_error_2d(exact, coefs, X, Y, Ni, LG, EQoLG;
    gauss_n = 5,
)
    local N_e = foldl(*, Ni)

    local ws, ps = gauss_quadrature_table[gauss_n]

    local coefs_ext = cat(coefs, 0.0, dims=1)

    local acc = 0.0

    for e in 1:N_e
        local LGe = LG[:, e]
        local Xe = X[LGe]
        local Ye = Y[LGe]

        local x2xis = x2xis_f(ps, Xe, Ye)
        local dx2xis = dx2xis_f(ps, Xe, Ye)

        for g_i = 1:gauss_n
            local p1 = ps[g_i]
            for g_j in 1:gauss_n
                local p2 = ps[g_j]
                local _x = x2xis[g_i][g_j]
                local _J = dx2xis[g_i][g_j]
                local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
                local diff = exact(x2xis[g_i][g_j])
                for i in 1:4
                    diff -= coefs_ext[EQoLG[i, e]] * phi[i]([p1, p2])
                end
                acc += J * ws[g_i] * ws[g_j] * diff * diff
            end
        end
    end

    sqrt(acc)
end

function instantiate_solution(coefs, phis, hi, Ni, EQoLG)
    coefs_ext = cat(coefs, 0.0, dims=1)
    call2 = (f, x, y) -> f(x, y)
    ( x -> begin
        ij = @. min(Int(floor(x / hi)), Ni-1)
        e = ij[2]*Ni[1] + ij[1] + 1
        pe = ij .* hi
        xii = @. (2 * (x - pe) / hi) - 1
        display(call2.(phis, xii...))
        foldl(+, (@. coefs_ext[EQoLG[:, e]] * call2(phis, xii...)))
    end
    )
end

end # module FiniteElements
