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
export phis_f, phi_derivs_f, x2xis_f, dx2xis_f
export build_small_vec_2d, build_vec_2d
export build_small_mat_2d, build_mat_2d

const phi = [
    (xi :: Vector{Float64} -> ((1 - xi[1])*(1 - xi[2]) / 4) :: Float64),
    (xi :: Vector{Float64} -> ((1 + xi[1])*(1 - xi[2]) / 4) :: Float64),
    (xi :: Vector{Float64} -> ((1 + xi[1])*(1 + xi[2]) / 4) :: Float64),
    (xi :: Vector{Float64} -> ((1 - xi[1])*(1 + xi[2]) / 4) :: Float64),
]

const phi_deriv = [
    ((xi...) -> ((- 0.25) * (1 - xi[2])) :: Float64)
    ((xi...) -> ((  0.25) * (1 - xi[2])) :: Float64)
    ((xi...) -> ((  0.25) * (1 + xi[2])) :: Float64)
    ((xi...) -> ((- 0.25) * (1 + xi[2])) :: Float64)
;;
    ((xi...) -> ((- 0.25) * (1 - xi[1])) :: Float64)
    ((xi...) -> ((- 0.25) * (1 + xi[1])) :: Float64)
    ((xi...) -> ((+ 0.25) * (1 + xi[1])) :: Float64)
    ((xi...) -> ((+ 0.25) * (1 - xi[1])) :: Float64)
]

const sdim = 2
const app = (f, xs...) -> f(xs...)

const phis_f = (ps :: AbstractVector{Float64}) ->
    map.(phi, ([
        [p1, p2]
        for p1 in ps, p2 in ps
    ],)) :: Vector{Matrix{Float64}}

# TODO: use array comprehension
const phi_derivs_f = (ps :: AbstractVector{Float64}) -> begin
    local phi_derivs = fill(0.0, (length(ps), length(ps), sdim, 4))

    for g_i in 1:length(ps)
        p_i = ps[g_i]
        for g_j in 1:length(ps)
            p_j = ps[g_j]
            for i in 1:sdim
                for j in 1:4
                    phi_derivs[g_i, g_j, i, j] = phi_deriv[j, i](p_i, p_j)
                end
            end
        end
    end
    phi_derivs
end :: Array{Float64, 4}

# TODO: Remove allocation (Vector creation)
const x2xis_f = (
    phis :: Vector{Matrix{Float64}},
    Xe :: AbstractVector{Float64}, Ye :: AbstractVector{Float64}
) -> begin
    local dot(xs, ys) = sum(xs .* ys)
    ((x, y) -> [x, y]).(dot(Xe, phis), dot(Ye, phis))
end :: Matrix{Vector{Float64}}

# TODO: use array comprehension
const dx2xis_f = (
    phi_derivs :: Array{Float64, 4},
    Xe :: AbstractVector{Float64}, Ye :: AbstractVector{Float64}
) -> begin
    dx2xis = fill(0.0, (size(phi_derivs)[1:2]..., length(Xe)))

    for p_i in 1:(size(phi_derivs)[1])
        for p_j in 1:(size(phi_derivs)[2])
            phi_derivs_1 = phi_derivs[p_i, p_j, 1, :]
            phi_derivs_2 = phi_derivs[p_i, p_j, 2, :]

            dx2xis[p_i, p_j, 1] = dot(Xe, phi_derivs_1)
            dx2xis[p_i, p_j, 2] = dot(Xe, phi_derivs_2)
            dx2xis[p_i, p_j, 3] = dot(Ye, phi_derivs_1)
            dx2xis[p_i, p_j, 4] = dot(Ye, phi_derivs_2)
        end
    end

    dx2xis
end :: Array{Float64, 3}

function build_small_mat_2d(
    alpha :: Float64, beta :: Float64,
    dx2xis :: Array{Float64, 3},
    ws :: Vector{Float64}, ps :: Vector{Float64}, gauss_n :: Int64,
) :: Matrix{Float64}
    dim = 4

    K_alpha1 = fill(0.0, (dim,dim))
    for i in 1:dim
        for j in 1:dim
            for g_i in 1:gauss_n
                p1 = ps[g_i]
                for g_j in 1:gauss_n
                    p2 = ps[g_j]
                    _J1 = dx2xis[g_i,g_j,1]
                    _J2 = dx2xis[g_i,g_j,2]
                    _J3 = dx2xis[g_i,g_j,3]
                    _J4 = dx2xis[g_i,g_j,4]
                    J = 1.0 / ((_J1 * _J4) - (_J2 * _J3))
                    H_1 = (_J4*_J4 + _J2*_J2)
                    H_2 = - (_J4*_J3 + _J2*_J1)
                    K_alpha1[i,j] += J * ws[g_i] * ws[g_j] * (
                        phi_deriv[j,1](p1, p2) * (
                            (H_1 * phi_deriv[i,1](p1, p2))
                            + (H_2 * phi_deriv[i,2](p1, p2))
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
                    _J1 = dx2xis[g_i,g_j,1]
                    _J2 = dx2xis[g_i,g_j,2]
                    _J3 = dx2xis[g_i,g_j,3]
                    _J4 = dx2xis[g_i,g_j,4]
                    J = 1.0 / ((_J1 * _J4) - (_J2 * _J3))
                    H_1 = - (_J4*_J3 + _J2*_J1)
                    H_2 = (_J3*_J3 + _J1*_J1)
                    K_alpha2[i,j] += J * ws[g_i] * ws[g_j] * (
                        phi_deriv[j,2](p1, p2) * (
                            (H_1 * phi_deriv[i,1](p1, p2))
                            + (H_2 * phi_deriv[i,2](p1, p2))
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
                    _J1 = dx2xis[g_i,g_j,1]
                    _J2 = dx2xis[g_i,g_j,2]
                    _J3 = dx2xis[g_i,g_j,3]
                    _J4 = dx2xis[g_i,g_j,4]
                    J = (_J1 * _J4) - (_J2 * _J3)
                    K_beta[i,j] += J * ws[g_i] * ws[g_j] *
                        (phi[j]([p1, p2])*phi[i]([p1, p2]))
                end
            end
        end
    end
    K_beta *= beta

    K_alpha1 + K_alpha2 + K_beta
end

function build_mat_2d(alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phi_derivs :: Array{Float64, 4},
    ws :: Vector{Float64}, ps :: Vector{Float64}, gauss_n :: Int64
) :: Matrix{Float64}

    local K = spzeros((m+1, m+1))
    for e in 1:N_e
        local LGe = LG[:, e]
        local Xe = X[LGe]
        local Ye = Y[LGe]

        local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

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

function build_small_vec_2d(f,
    x2xis :: Matrix{Vector{Float64}},
    dx2xis :: Array{Float64, 3},
    phis :: Vector{Matrix{Float64}},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}
    local dim = 4

    local F = fill(0.0, (dim,))
    for i in 1:dim
        for g_i in 1:gauss_n
            for g_j in 1:gauss_n
                local _x = x2xis[g_i,g_j]
                local _J1 = dx2xis[g_i,g_j,1]
                local _J2 = dx2xis[g_i,g_j,2]
                local _J3 = dx2xis[g_i,g_j,3]
                local _J4 = dx2xis[g_i,g_j,4]
                local J = (_J1 * _J4) - (_J2 * _J3)
                F[i] += J * ws[g_i] * ws[g_j] * (
                    f(_x)*phis[i][g_i,g_j]
                )
            end
        end
    end

    F
end

function build_vec_2d(f,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Vector{Matrix{Float64}},
    phi_derivs :: Array{Float64, 4},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}

    local F = fill(0.0, (m+1,))
    for e in 1:N_e
        local LGe = LG[:, e]
        local Xe = X[LGe]
        local Ye = Y[LGe]

        local x2xis = x2xis_f(phis, Xe, Ye)
        local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

        local F_e = build_small_vec_2d(f, x2xis,
            dx2xis, phis, ws, gauss_n
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

function build_LG(Ni :: AbstractVector{Int64}) :: AbstractMatrix{Int64}
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

function build_EQ(Ni :: AbstractVector{Int64}) :: Tuple{Int64, Vector{Int64}}
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

function finite_elements_setup(ex :: Example,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    Ni :: AbstractVector{Int64}
;
    gauss_n :: Int64 = 5
)
    finite_elements_setup(ex.f, ex.alpha, ex.beta, X, Y, Ni,
        gauss_n=gauss_n
    )
end

function finite_elements_setup(f,
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    Ni :: AbstractVector{Int64}
;
    gauss_n :: Int64 = 5
)
    local N_e = foldl(*, Ni)

    local LG = build_LG(Ni)

    local m, EQ = build_EQ(Ni)

    local EQoLG = EQ[LG]

    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)

    local K = build_mat_2d(alpha, beta, X, Y, N_e, LG, EQoLG, m,
        phi_derivs,
        ws, ps, gauss_n,
    )
    local F = build_vec_2d(f, X, Y, N_e, LG, EQoLG, m,
        phis, phi_derivs,
        ws, gauss_n,
    )

    K, F, LG, EQoLG, m
end

function generate_space(
    hi :: AbstractVector{Float64}, Ni :: AbstractVector{Int64}
;
    noise :: Bool = false
) :: Tuple{AbstractVector{Float64}, AbstractVector{Float64}}
    X = repeat(0.0:(hi[1]):1.0, Ni[2]+1)
    Y = cat(
        ((x->repeat(x:x, Ni[1]+1)).(0.0:(hi[2]):1.0))...,
        dims=1
    )

    if noise
        scale = hi ./ 4
        for i in 2:Ni[2]
            for j in 2:Ni[1]
                idx = (i-1) * (Ni[1]+1) + (j-1) + 1
                X[idx] += scale[1] * 2*(rand(Float64) - 0.5)
                Y[idx] += scale[2] * 2*(rand(Float64) - 0.5)
            end
        end
    end

    X, Y
end

function gauss_error_2d(exact,
    coefs,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    Ni :: AbstractVector{Int64}, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}
;
    gauss_n :: Int64 = 5,
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
