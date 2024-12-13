module FiniteElements

include("../examples2d.jl")
using .Examples2d: Example, example

include("../common.jl")
using .Common: gauss_quadrature_table

using LinearAlgebra: dot

using SparseArrays: spzeros

import Random

export finite_elements_setup
export gauss_quadrature_table
export generate_space
export gauss_error_2d
export phi, phi_deriv, instantiate_solution
export Example, example

const phi = [
    ((xi...) -> ((1 - xi[1])*(1 - xi[2]) / 4) :: Float64),
    ((xi...) -> ((1 + xi[1])*(1 - xi[2]) / 4) :: Float64),
    ((xi...) -> ((1 + xi[1])*(1 + xi[2]) / 4) :: Float64),
    ((xi...) -> ((1 - xi[1])*(1 + xi[2]) / 4) :: Float64),
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

const dim = 4
const sdim = 2
const app = (f, xs...) -> f(xs...)

const phis_f = (ps :: AbstractVector{Float64}) -> begin
    local phis = fill(0.0, (length(ps), length(ps), length(phi)))
    for (g_i, p1) in enumerate(ps),
        (g_j, p2) in enumerate(ps),
        (i, f_phi) in enumerate(phi)
        phis[g_i,g_j,i] = f_phi(p1, p2)
    end
    phis
end :: Array{Float64, 3}

const phi_derivs_f = (ps :: AbstractVector{Float64}) -> begin
    phi_derivs = fill(0.0, (length(ps), sdim, dim))
    for g_i in 1:length(ps),
        j in 1:sdim,
        i in 1:dim,
        phi_derivs[g_i,j,i] = phi_deriv[i,j](ps[g_i], ps[g_i])
    end
    phi_derivs
end :: Array{Float64, 3}

const x2xis_f = (
    phis :: Array{Float64, 3},
    Xe :: AbstractVector{Float64}, Ye :: AbstractVector{Float64}
) -> begin
    local x2xis = fill(0.0, (size(phis)[1:2]..., sdim))
    for i in 1:(size(phis)[1]),
        j in 1:(size(phis)[2])
        local phis_ = view(phis, i, j, :)
        x2xis[i,j,1] = dot(Xe, phis_)
        x2xis[i,j,2] = dot(Ye, phis_)
    end
    x2xis
end :: Array{Float64, 3}

const x2xis_f_ref! = (
    ref_out :: Ref{Array{Float64, 3}},
    phis :: Array{Float64, 3},
    Xe :: AbstractVector{Float64}, Ye :: AbstractVector{Float64}
) -> begin
    for i in 1:(size(phis)[1]),
        j in 1:(size(phis)[2])
        local phis_ = view(phis, i, j, :)
        ref_out[][i,j,1] = dot(Xe, phis_)
        ref_out[][i,j,2] = dot(Ye, phis_)
    end
    ref_out
end :: Ref{Array{Float64, 3}}

const dx2xis_f = (
    phi_derivs :: Array{Float64, 3},
    Xe :: AbstractVector{Float64}, Ye :: AbstractVector{Float64}
) -> begin
    dx2xis = fill(0.0, (size(phi_derivs)[1], size(phi_derivs)[1], length(Xe)))

    for p_i in 1:(size(phi_derivs)[1])
        phi_derivs_i = view(phi_derivs, p_i, 2, :)
        for p_j in 1:(size(phi_derivs)[1])
            phi_derivs_j = view(phi_derivs, p_j, 1, :)

            dx2xis[p_i, p_j, 1] = dot(Xe, phi_derivs_j)
            dx2xis[p_i, p_j, 2] = dot(Xe, phi_derivs_i)
            dx2xis[p_i, p_j, 3] = dot(Ye, phi_derivs_j)
            dx2xis[p_i, p_j, 4] = dot(Ye, phi_derivs_i)
        end
    end

    dx2xis
end :: Array{Float64, 3}

const dx2xis_f_ref! = (
    ref_out :: Ref{Array{Float64, 3}},
    phi_derivs :: Array{Float64, 3},
    Xe :: AbstractVector{Float64}, Ye :: AbstractVector{Float64}
) -> begin
    for p_i in 1:(size(phi_derivs)[1])
        phi_derivs_i = view(phi_derivs, p_i, 2, :)
        for p_j in 1:(size(phi_derivs)[1])
            phi_derivs_j = view(phi_derivs, p_j, 1, :)

            ref_out[][p_i, p_j, 1] = dot(Xe, phi_derivs_j)
            ref_out[][p_i, p_j, 2] = dot(Xe, phi_derivs_i)
            ref_out[][p_i, p_j, 3] = dot(Ye, phi_derivs_j)
            ref_out[][p_i, p_j, 4] = dot(Ye, phi_derivs_i)
        end
    end
    ref_out
end :: Ref{Array{Float64, 3}}

function build_small_vec_2d(
    f :: Function,
    x2xis :: Array{Float64, 3},
    dx2xis :: Array{Float64, 3},
    phis :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}

    local F = fill(0.0, (dim,))
    for g_i in 1:gauss_n
        for g_j in 1:gauss_n
            local _x = view(x2xis, g_i, g_j, :)
            local _J = view(dx2xis, g_i, g_j, :)
            local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
            local pre_calc = J * ws[g_i] * ws[g_j] * f(_x[1], _x[2])
            local phis_ = view(phis, g_i, g_j, :)
            for i in 1:dim
                F[i] += pre_calc * phis_[i]
            end
        end
    end
    F
end

function build_small_vec_2d_ref!(
    ref_out :: Ref{Vector{Float64}},
    f :: Function,
    x2xis :: Array{Float64, 3},
    dx2xis :: Array{Float64, 3},
    phis :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Ref{Vector{Float64}}

    fill!(ref_out[], 0.0)
    for g_i in 1:gauss_n
        for g_j in 1:gauss_n
            local _x = view(x2xis, g_i, g_j, :)
            local _J = view(dx2xis, g_i, g_j, :)
            local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
            local pre_calc = J * ws[g_i] * ws[g_j] * f(_x[1], _x[2])
            local phis_ = view(phis, g_i, g_j, :)
            for i in 1:dim
                ref_out[][i] += pre_calc * phis_[i]
            end
        end
    end
    ref_out
end

function build_vec_2d(
    f :: Function,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}

    local F = fill(0.0, (m+1,))
    for e in 1:N_e
        local LGe = view(LG, :, e)
        local Xe = view(X, LGe)
        local Ye = view(Y, LGe)

        local x2xis = x2xis_f(phis, Xe, Ye)
        local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

        local F_e = build_small_vec_2d(
            f,
            x2xis, dx2xis,
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += F_e[i]
        end
    end
    F[begin:end-1]
end

function build_vec_2d_ref(
    f :: Function,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}

    local F = fill(0.0, (m+1,))
    local ref_Fe = Ref(fill(0.0, (dim,)))
    local ref_x2xis = Ref(fill(0.0, (size(phis)[1:2]..., sdim)))
    local ref_dx2xis = Ref(fill(0.0, (size(phi_derivs)[1], size(phi_derivs)[1], dim)))
    for e in 1:N_e
        local LGe = view(LG, :, e)
        local Xe = view(X, LGe)
        local Ye = view(Y, LGe)

        x2xis_f_ref!(ref_x2xis, phis, Xe, Ye)
        dx2xis_f_ref!(ref_dx2xis, phi_derivs, Xe, Ye)

        build_small_vec_2d_ref!(
            ref_Fe,
            f,
            ref_x2xis[], ref_dx2xis[],
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += ref_Fe[][i]
        end
    end
    F[begin:end-1]
end

function build_vec_2d_iter(
    f :: Function,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}

    local iter = FiniteElementIter(
        X, Y,
        N_e, LG,
        phis, phi_derivs,
    )

    local F = fill(0.0, (m+1,))
    for (e, x2xis, dx2xis) in iter
        local F_e = build_small_vec_2d(
            f,
            x2xis, dx2xis,
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += F_e[i]
        end
    end
    F[begin:end-1]
end

function build_vec_2d_iterref(
    f :: Function,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Vector{Float64}

    local iter = FiniteElementIterRef(
        X, Y,
        N_e, LG,
        phis, phi_derivs,
    )

    local F = fill(0.0, (m+1,))
    for (e, ref_x2xis, ref_dx2xis, ref_Ke, ref_Fe) in iter
        local x2xis = ref_x2xis[]
        local dx2xis = ref_dx2xis[]

        build_small_vec_2d_ref!(
            ref_Fe,
            f,
            x2xis, dx2xis,
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += ref_Fe[][i]
        end
    end
    F[begin:end-1]
end

function build_small_mat_2d(
    alpha :: Float64, beta :: Float64,
    dx2xis :: Array{Float64, 3},
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64,
) :: Matrix{Float64}

    K = fill(0.0, (dim,dim))

    for g_i in 1:gauss_n
        local phi_derivs_i = view(phi_derivs, g_i, 2, :)
        for g_j in 1:gauss_n
            local phis_ = view(phis, g_i, g_j, :)

            local _J = view(dx2xis, g_i, g_j, :)
            local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
            local invJ = 1.0 / J

            local H1_1 = (_J[4]*_J[4] + _J[2]*_J[2])
            local H1_2 = - (_J[4]*_J[3] + _J[2]*_J[1])

            local H2_1 = - (_J[4]*_J[3] + _J[2]*_J[1])
            local H2_2 = (_J[3]*_J[3] + _J[1]*_J[1])

            local alpha_pre_calc = alpha * invJ * ws[g_i] * ws[g_j]
            local beta_pre_calc = beta * J * ws[g_i] * ws[g_j]
            local phi_derivs_j = view(phi_derivs, g_j, 1, :)
            for i in 1:dim
                local alpha1_pre_calc = alpha_pre_calc * (
                    (H1_1 * phi_derivs_j[i])
                    + (H1_2 * phi_derivs_i[i])
                )
                local alpha2_pre_calc = alpha_pre_calc * (
                    (H2_1 * phi_derivs_j[i])
                    + (H2_2 * phi_derivs_i[i])
                )
                for j in 1:dim
                    local kij = (
                        (alpha1_pre_calc * phi_derivs_j[j])
                        + (alpha2_pre_calc * phi_derivs_i[j])
                        + (beta_pre_calc * (phis_[j]*phis_[i]))
                    )
                    K[i,j] += kij
                end
            end
        end
    end
    K
end

function build_small_mat_2d_ref!(
    ref_out :: Ref{Matrix{Float64}},
    alpha :: Float64, beta :: Float64,
    dx2xis :: Array{Float64, 3},
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64,
) :: Ref{Matrix{Float64}}

    fill!(ref_out[], 0.0)

    for g_i in 1:gauss_n
        local phi_derivs_i = view(phi_derivs, g_i, 2, :)
        for g_j in 1:gauss_n
            local phis_ = view(phis, g_i, g_j, :)

            local _J = view(dx2xis, g_i, g_j, :)
            local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
            local invJ = 1.0 / J

            local H1_1 = (_J[4]*_J[4] + _J[2]*_J[2])
            local H1_2 = - (_J[4]*_J[3] + _J[2]*_J[1])

            local H2_1 = - (_J[4]*_J[3] + _J[2]*_J[1])
            local H2_2 = (_J[3]*_J[3] + _J[1]*_J[1])

            local alpha_pre_calc = alpha * invJ * ws[g_i] * ws[g_j]
            local beta_pre_calc = beta * J * ws[g_i] * ws[g_j]
            local phi_derivs_j = view(phi_derivs, g_j, 1, :)
            for i in 1:dim
                local alpha1_pre_calc = alpha_pre_calc * (
                    (H1_1 * phi_derivs_j[i])
                    + (H1_2 * phi_derivs_i[i])
                )
                local alpha2_pre_calc = alpha_pre_calc * (
                    (H2_1 * phi_derivs_j[i])
                    + (H2_2 * phi_derivs_i[i])
                )
                for j in 1:dim
                    local kij = (
                        (alpha1_pre_calc * phi_derivs_j[j])
                        + (alpha2_pre_calc * phi_derivs_i[j])
                        + (beta_pre_calc * (phis_[j]*phis_[i]))
                    )
                    ref_out[][i,j] += kij
                end
            end
        end
    end
    ref_out
end

function build_mat_2d(
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: AbstractMatrix{Float64}

    local K = spzeros((m+1, m+1))
    for e in 1:N_e
        local LGe = view(LG, :, e)
        local Xe = view(X, LGe)
        local Ye = view(Y, LGe)

        local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

        local K_e = build_small_mat_2d(
            alpha, beta,
            dx2xis,
            phis, phi_derivs,
            ws, gauss_n,
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += K_e[i, j]
            end
        end
    end
    K[begin:end-1, begin:end-1]
end

function build_mat_2d_ref(
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: AbstractMatrix{Float64}

    local K = spzeros((m+1, m+1))
    local ref_Ke = Ref(fill(0.0, (dim, dim)))
    local ref_dx2xis = Ref(fill(0.0, (size(phi_derivs)[1], size(phi_derivs)[1], dim)))
    for e in 1:N_e
        local LGe = view(LG, :, e)
        local Xe = view(X, LGe)
        local Ye = view(Y, LGe)

        dx2xis_f_ref!(ref_dx2xis, phi_derivs, Xe, Ye)

        build_small_mat_2d_ref!(
            ref_Ke,
            alpha, beta,
            ref_dx2xis[],
            phis, phi_derivs,
            ws, gauss_n,
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += ref_Ke[][i, j]
            end
        end
    end
    K[begin:end-1, begin:end-1]
end

function build_mat_2d_iter(
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: AbstractMatrix{Float64}

    local iter = FiniteElementIter(
        X, Y,
        N_e, LG,
        phis, phi_derivs,
    )

    local K = spzeros((m+1, m+1))
    for (e, x2xis, dx2xis) in iter
        local K_e = build_small_mat_2d(
            alpha, beta,
            dx2xis,
            phis, phi_derivs,
            ws, gauss_n,
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += K_e[i, j]
            end
        end
    end
    K[begin:end-1, begin:end-1]
end

function build_mat_2d_iterref(
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: AbstractMatrix{Float64}

    local iter = FiniteElementIterRef(
        X, Y,
        N_e, LG,
        phis, phi_derivs,
    )

    local K = spzeros((m+1, m+1))
    for (e, ref_x2xis, ref_dx2xis, ref_Ke, ref_Fe) in iter
        local x2xis = ref_x2xis[]
        local dx2xis = ref_dx2xis[]
        build_small_mat_2d_ref!(
            ref_Ke,
            alpha, beta,
            dx2xis,
            phis, phi_derivs,
            ws, gauss_n,
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += ref_Ke[][i, j]
            end
        end
    end
    K[begin:end-1, begin:end-1]
end

function build_vec_mat_2d(
    f :: Function,
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Tuple{Vector{Float64}, AbstractMatrix{Float64}}

    local K = spzeros((m+1, m+1))
    local F = fill(0.0, (m+1,))
    for e in 1:N_e
        local LGe = view(LG, :, e)
        local Xe = view(X, LGe)
        local Ye = view(Y, LGe)

        local x2xis = x2xis_f(phis, Xe, Ye)
        local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

        local K_e = build_small_mat_2d(
            alpha, beta,
            dx2xis,
            phis, phi_derivs,
            ws, gauss_n,
        )
        local F_e = build_small_vec_2d(
            f,
            x2xis, dx2xis,
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += F_e[i]
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += K_e[i, j]
            end
        end
    end

    (F[begin:end-1], K[begin:end-1,begin:end-1])
end

function build_vec_mat_2d_ref(
    f :: Function,
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Tuple{Vector{Float64}, AbstractMatrix{Float64}}

    local K = spzeros((m+1, m+1))
    local F = fill(0.0, (m+1,))
    local ref_x2xis = Ref(fill(0.0, (size(phis)[1:2]..., sdim)))
    local ref_dx2xis = Ref(fill(0.0, (size(phi_derivs)[1], size(phi_derivs)[1], dim)))
    for e in 1:N_e
        local LGe = view(LG, :, e)
        local Xe = view(X, LGe)
        local Ye = view(Y, LGe)

        x2xis_f_ref!(ref_x2xis, phis, Xe, Ye)
        dx2xis_f_ref!(ref_dx2xis, phi_derivs, Xe, Ye)

        local K_e = build_small_mat_2d(
            alpha, beta,
            ref_dx2xis[],
            phis, phi_derivs,
            ws, gauss_n,
        )
        local F_e = build_small_vec_2d(
            f,
            ref_x2xis[], ref_dx2xis[],
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += F_e[i]
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += K_e[i, j]
            end
        end
    end

    (F[begin:end-1], K[begin:end-1,begin:end-1])
end

function build_vec_mat_2d_iter(
    f :: Function,
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Tuple{Vector{Float64}, AbstractMatrix{Float64}}

    local iter = FiniteElementIter(
        X, Y,
        N_e, LG,
        phis, phi_derivs,
    )

    local K = spzeros((m+1, m+1))
    local F = fill(0.0, (m+1,))
    for (e, x2xis, dx2xis) in iter
        local K_e = build_small_mat_2d(
            alpha, beta,
            dx2xis,
            phis, phi_derivs,
            ws, gauss_n,
        )
        local F_e = build_small_vec_2d(
            f,
            x2xis, dx2xis,
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += F_e[i]
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += K_e[i, j]
            end
        end
    end

    (F[begin:end-1], K[begin:end-1,begin:end-1])
end

function build_vec_mat_2d_iterref(
    f :: Function,
    alpha :: Float64, beta :: Float64,
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    N_e :: Int64, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}, m :: Int64,
    phis :: Array{Float64, 3},
    phi_derivs :: Array{Float64, 3},
    ws :: Vector{Float64}, gauss_n :: Int64
) :: Tuple{Vector{Float64}, AbstractMatrix{Float64}}

    local iter = FiniteElementIterRef(
        X, Y,
        N_e, LG,
        phis, phi_derivs,
    )

    local K = spzeros((m+1, m+1))
    local F = fill(0.0, (m+1,))
    for (e, ref_x2xis, ref_dx2xis, ref_Ke, ref_Fe) in iter
        local x2xis = ref_x2xis[]
        local dx2xis = ref_dx2xis[]
        build_small_mat_2d_ref!(
            ref_Ke,
            alpha, beta,
            dx2xis,
            phis, phi_derivs,
            ws, gauss_n,
        )
        build_small_vec_2d_ref!(
            ref_Fe,
            f,
            x2xis, dx2xis,
            phis, ws, gauss_n
        )

        local EQoLG_ = view(EQoLG, :, e)
        for i in 1:dim
            F[EQoLG_[i]] += ref_Fe[][i]
            for j in 1:dim
                K[EQoLG_[i], EQoLG_[j]] += ref_Ke[][i, j]
            end
        end
    end

    (F[begin:end-1], K[begin:end-1,begin:end-1])
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

    local F, K = build_vec_mat_2d_ref(
        f, alpha, beta,
        X, Y,
        N_e, LG, EQoLG, m,
        phis, phi_derivs,
        ws, gauss_n,
    )

    K, F, LG, EQoLG, m
end

function generate_space(
    hi :: AbstractVector{Float64}, Ni :: AbstractVector{Int64}
;
    noise :: Bool = false,
    seed :: Union{Int64, Nothing} = nothing,
) :: Tuple{AbstractVector{Float64}, AbstractVector{Float64}}
    X = repeat(0.0:(hi[1]):1.0, Ni[2]+1)
    Y = cat(
        ((x->repeat(x:x, Ni[1]+1)).(0.0:(hi[2]):1.0))...,
        dims=1
    )

    if noise
        if seed !== nothing
            Random.seed!(seed)
        end
        scale = hi ./ 4
        for i in 2:Ni[2]
            for j in 2:Ni[1]
                idx = (i-1) * (Ni[1]+1) + (j-1) + 1
                X[idx] += scale[1] * 2*(Random.rand(Float64) - 0.5)
                Y[idx] += scale[2] * 2*(Random.rand(Float64) - 0.5)
            end
        end
    end

    X, Y
end

function gauss_error_2d(
    exact :: Function,
    coefs :: AbstractVector{Float64},
    X :: AbstractVector{Float64}, Y :: AbstractVector{Float64},
    Ni :: AbstractVector{Int64}, LG :: AbstractMatrix{Int64},
    EQoLG :: Matrix{Int64}
;
    gauss_n :: Int64 = 5,
)
    local N_e = foldl(*, Ni)

    local ws, ps = gauss_quadrature_table[gauss_n]

    local phis = phis_f(ps)
    local phi_derivs = phi_derivs_f(ps)

    local coefs_ext = cat(coefs, 0.0, dims=1)

    local acc = 0.0

    for e in 1:N_e
        local LGe = LG[:, e]
        local Xe = X[LGe]
        local Ye = Y[LGe]

        local x2xis = x2xis_f(phis, Xe, Ye)
        local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

        for g_i = 1:gauss_n
            for g_j in 1:gauss_n
                local _x = view(x2xis, g_i, g_j, :)
                local _J = view(dx2xis, g_i, g_j, :)
                local J = (_J[1] * _J[4]) - (_J[2] * _J[3])
                local diff = exact(_x...)
                for i in 1:4
                    diff -= coefs_ext[EQoLG[i, e]] * phis[g_i, g_j, i]
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

const FiniteElementState = Int64
const FiniteElementItem = Tuple{
    Int64, Array{Float64, 3}, Array{Float64, 3}
}

struct FiniteElementIter
    X :: AbstractVector{Float64}
    Y :: AbstractVector{Float64}
    N_e :: Int64
    LG :: AbstractMatrix{Int64}
    phis :: Array{Float64, 3}
    phi_derivs :: Array{Float64, 3}

    FiniteElementIter(
        X :: AbstractVector{Float64},
        Y :: AbstractVector{Float64},
        N_e :: Int64,
        LG :: AbstractMatrix{Int64},
        phis :: Array{Float64, 3},
        phi_derivs :: Array{Float64, 3},
    ) = begin
        new(
            X, Y,
            N_e, LG,
            phis, phi_derivs,
        )
    end
end

const Base.iterate(
    iter :: FiniteElementIter,
    state :: Int = 1,
) = (!(state <= length(iter))) ? nothing : begin
    local next_state = state + 1
    local item = getindex(iter, state)
    (item, next_state)
end :: Union{Tuple{FiniteElementItem, FiniteElementState}, Nothing}

const Base.IteratorSize(::Type{FiniteElementIter}) = Base.HasLength()
const Base.length(iter :: FiniteElementIter) = iter.N_e
const Base.IteratorEltype(::Type{FiniteElementIter}) = Base.HasEltype()
const Base.eltype(::Type{FiniteElementIter}) = FiniteElementItem
const Base.isdone(iter :: FiniteElementIter, state :: Int64 = 0) = Base.length(iter) == state

const Base.getindex(iter :: FiniteElementIter, i :: FiniteElementState) = begin
    local LGe = view(iter.LG, :, i)
    local Xe = view(iter.X, LGe)
    local Ye = view(iter.Y, LGe)
    local x2xis  = x2xis_f(iter.phis, Xe, Ye)
    local dx2xis = dx2xis_f(iter.phi_derivs, Xe, Ye)
    (i, x2xis, dx2xis)
end :: FiniteElementItem
const Base.firstindex(iter :: FiniteElementIter) = 1
const Base.lastindex(iter :: FiniteElementIter) = length(iter)
const Base.IndexStyle(::Type{FiniteElementIter}) = Base.IndexLinear()

const FiniteElementItemRef = Tuple{
    Int64, Ref{Array{Float64, 3}}, Ref{Array{Float64, 3}},
    Ref{Matrix{Float64}}, Ref{Vector{Float64}},
}

struct FiniteElementIterRef
    iter :: FiniteElementIter
    x2xis :: Ref{Array{Float64, 3}}
    dx2xis :: Ref{Array{Float64, 3}}
    Ke :: Ref{Matrix{Float64}}
    Fe :: Ref{Vector{Float64}}

    FiniteElementIterRef(
        iter :: FiniteElementIter,
    ) = begin
        local x2xis = Ref(fill(0.0, (size(iter.phis)[1:2]..., sdim)))
        local dx2xis = Ref(fill(0.0, (size(iter.phi_derivs)[1], size(iter.phi_derivs)[1], dim)))
        local Ke = Ref(fill(0.0, (dim, dim)))
        local Fe = Ref(fill(0.0, (dim,)))
        new(iter, x2xis, dx2xis, Ke, Fe)
    end
    FiniteElementIterRef(
        X :: AbstractVector{Float64},
        Y :: AbstractVector{Float64},
        N_e :: Int64,
        LG :: AbstractMatrix{Int64},
        phis :: Array{Float64, 3},
        phi_derivs :: Array{Float64, 3},
    ) = FiniteElementIterRef(
        FiniteElementIter(
            X, Y,
            N_e, LG,
            phis, phi_derivs,
        )
    )
end

function getindex!(self :: FiniteElementIterRef, e :: FiniteElementState) :: FiniteElementItemRef
    local iter = self.iter
    local LGe = view(iter.LG, :, e)
    local Xe = view(iter.X, LGe)
    local Ye = view(iter.Y, LGe)
    x2xis_f_ref!(self.x2xis, iter.phis, Xe, Ye)
    dx2xis_f_ref!(self.dx2xis, iter.phi_derivs, Xe, Ye)
    (e, self.x2xis, self.dx2xis, self.Ke, self.Fe)
end

const Base.iterate(
    iter :: FiniteElementIterRef,
    state :: Int = 1,
) = (!(state <= length(iter))) ? nothing : begin
    local next_state = state + 1
    local item = getindex!(iter, state)
    (item, next_state)
end :: Union{Tuple{FiniteElementItemRef, FiniteElementState}, Nothing}

const Base.IteratorSize(::Type{FiniteElementIterRef}) = Base.IteratorSize(Type{FiniteElementIter})
const Base.length(iter :: FiniteElementIterRef) = Base.length(iter.iter)
const Base.IteratorEltype(::Type{FiniteElementIterRef}) = Base.IteratorEltype(Type{FiniteElementIter})
const Base.eltype(::Type{FiniteElementIterRef}) = Base.eltype(Type{FiniteElementIter})
const Base.isdone(iter :: FiniteElementIterRef, state :: Int64 = 0) = Base.isdone(iter.iter, state)

end # module FiniteElements
