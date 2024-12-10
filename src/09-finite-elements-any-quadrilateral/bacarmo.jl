# Adapted from
# https://github.com/bacarmo/Elementos_Finitos/blob/dba5361b423a32933118b3fca65edb93fe1cc425/2024_02/Estacionario_2D_equacao1_Pluto.jl
module Bacarmo

using SparseArrays, LinearAlgebra

### Begin Hack
include("../common.jl")
using .Common: gauss_quadrature_table
function legendre(n)
    gauss_quadrature_table[n]
end
### End Hack

function ϕ(ξ₁::Float64, ξ₂::Float64) :: Vector{Float64}
    [(1-ξ₁)*(1-ξ₂)/4, (1+ξ₁)*(1-ξ₂)/4, (1+ξ₁)*(1+ξ₂)/4, (1-ξ₁)*(1+ξ₂)/4]
end

function ∂ϕ_∂ξ₁(ξ₂::Float64) :: Vector{Float64}
    [-(1-ξ₂)/4, (1-ξ₂)/4, (1+ξ₂)/4, -(1+ξ₂)/4]
end

function ∂ϕ_∂ξ₂(ξ₁::Float64) :: Vector{Float64}
    [-(1-ξ₁)/4, -(1+ξ₁)/4, (1+ξ₁)/4, (1-ξ₁)/4]
end

function monta_Fᵉ_quadrilatero!(Fᵉ::Vector{Float64}, f::Function,
                                X1e::Vector{Float64}, X2e::Vector{Float64},
                                P::Vector{Float64}, W::Vector{Float64})
    fill!(Fᵉ, 0.0)

    for i in 1:length(P)
        ξ₁ = P[i]

        vec_∂ϕ_∂ξ₂ = ∂ϕ_∂ξ₂(ξ₁)

        for j in 1:length(P)
            ξ₂ = P[j]

            vec_ϕ = ϕ(ξ₁,ξ₂)

            vec_∂ϕ_∂ξ₁ = ∂ϕ_∂ξ₁(ξ₂)

            x₁ = dot(X1e, vec_ϕ)
            x₂ = dot(X2e, vec_ϕ)

            detJ = dot(X1e, vec_∂ϕ_∂ξ₁) * dot(X2e, vec_∂ϕ_∂ξ₂) -
                   dot(X1e, vec_∂ϕ_∂ξ₂) * dot(X2e, vec_∂ϕ_∂ξ₁)
            @assert detJ > 0 "O determinante jacobiano deve ser positivo"

            for a in 1:4
                Fᵉ[a] += W[i] * W[j] * f(x₁, x₂) * vec_ϕ[a] * detJ
            end
        end
    end
end

function monta_Kᵉ_quadrilatero!(Kᵉ::Matrix{Float64}, α::Float64, β::Float64,
                                X1e::Vector{Float64}, X2e::Vector{Float64},
                                P::Vector{Float64}, W::Vector{Float64})
    fill!(Kᵉ, 0.0)

    for i in 1:length(P)
        ξ₁ = P[i]

        vec_∂ϕ_∂ξ₂ = ∂ϕ_∂ξ₂(ξ₁)

        for j in 1:length(P)
            ξ₂ = P[j]

            vec_ϕ = ϕ(ξ₁,ξ₂)

            vec_∂ϕ_∂ξ₁ = ∂ϕ_∂ξ₁(ξ₂)

            ∂x₁_∂ξ₁ = dot(X1e, vec_∂ϕ_∂ξ₁)
            ∂x₁_∂ξ₂ = dot(X1e, vec_∂ϕ_∂ξ₂)
            ∂x₂_∂ξ₁ = dot(X2e, vec_∂ϕ_∂ξ₁)
            ∂x₂_∂ξ₂ = dot(X2e, vec_∂ϕ_∂ξ₂)

            detJ = ∂x₁_∂ξ₁ * ∂x₂_∂ξ₂ - ∂x₁_∂ξ₂ * ∂x₂_∂ξ₁
            @assert detJ > 0 "O determinante jacobiano deve ser positivo"

            HᵀH₁₁ = ∂x₂_∂ξ₂^2 + ∂x₁_∂ξ₂^2
            HᵀH₁₂ = -∂x₁_∂ξ₁ * ∂x₁_∂ξ₂ - ∂x₂_∂ξ₁ * ∂x₂_∂ξ₂
            HᵀH₂₂ = ∂x₂_∂ξ₁^2 + ∂x₁_∂ξ₁^2

            for b in 1:4
                for a in 1:4
                    Kᵉ[a, b] += W[i] * W[j] * (
                        (α / detJ) * (
                            vec_∂ϕ_∂ξ₁[b] *
                            (HᵀH₁₁ * vec_∂ϕ_∂ξ₁[a] + HᵀH₁₂ * vec_∂ϕ_∂ξ₂[a]) +
                            vec_∂ϕ_∂ξ₂[b] *
                            (HᵀH₁₂ * vec_∂ϕ_∂ξ₁[a] + HᵀH₂₂ * vec_∂ϕ_∂ξ₂[a])
                        ) + β * vec_ϕ[b] * vec_ϕ[a] * detJ
                    )
                end
            end
        end
    end
end

function monta_F_quadrilatero(
    f::Function, X₁::AbstractArray{Float64}, X₂::AbstractArray{Float64},
    m::Int64, EQ::Vector{Int64}, LG::AbstractMatrix{Int64}) :: Vector{Float64}

    ne = size(LG,2)

    P, W = legendre(5)

    Fᵉ = zeros(4)
    F = zeros(m+1)

    for e in 1:ne
        idx = LG[:,e]

        X1e = X₁[idx]
        X2e = X₂[idx]

        idx = EQ[idx]

        monta_Fᵉ_quadrilatero!(Fᵉ, f, X1e, X2e, P, W)

        for a in 1:4
            F[idx[a]] += Fᵉ[a]
        end
    end

    return F[1:m]
end

function monta_K_quadrilatero(α::Float64, β::Float64,
                              X₁::AbstractArray{Float64}, X₂::AbstractArray{Float64},
                              m::Int64, EQ::Vector{Int64}, LG::AbstractMatrix{Int64}) ::
                              SparseMatrixCSC{Float64, Int64}
    ne = size(LG, 2)

    P, W = legendre(2)

    Kᵉ = zeros(4,4)

    K = spzeros(m+1, m+1)

    for e = 1:ne
        idx = LG[:,e]

        X1e = X₁[idx]
        X2e = X₂[idx]

        idx = EQ[idx]

        monta_Kᵉ_quadrilatero!(Kᵉ, α, β, X1e, X2e, P, W)

        for b = 1:4
            j = idx[b]
            for a = 1:4
                i = idx[a]
                K[i,j] += Kᵉ[a,b]
            end
        end
    end

    return K[1:m, 1:m]
end

end # Bacarmo
