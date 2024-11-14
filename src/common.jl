module Common

using LinearAlgebra: SymTridiagonal

export build_tridiagonal
function build_tridiagonal(main, side, dim)
    mains = fill(main, (dim,))
    sides = fill(side, (dim-1,))
    SymTridiagonal(mains, sides)
end

export n_points_from_to
function n_points_from_to(n; from=0, to=1, i_start=1, i_end=n)
    inter = (i_start:i_end) ./ (n+1)
    xs = (to - from) * inter .+ from
    xs
end

export calc_error
function calc_error(actual, expected)
    # norm(actual - expected) / norm(expected)
    maximum(abs.(actual - expected))
end

export gauss_quadrature_table
# Copied from: https://pomax.github.io/bezierinfo/legendre-gauss.html
const gauss_quadrature_table = [
    # 1
    ([2.0],
        [0.0]),
    # 2
    ([1.0, 1.0],
        [-0.5773502691896257, 0.5773502691896257]),
    # 3
    ([0.5555555555555556, 0.8888888888888888, 0.5555555555555556],
        [-0.7745966692414834, 0.0000000000000000, 0.7745966692414834]),
    # 4
    ([0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538],
        [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]),
    # 5
    ([0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891],
        [-0.9061798459386640, -0.5384693101056831, 0.0000000000000000, 0.5384693101056831, 0.9061798459386640]),
    # 6
    ([0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704],
        [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521]),
]

export gauss_quadrature
function gauss_quadrature(func, n)
    @assert n <= 6
    ws, xs = gauss_quadrature_table[n]

    fxs = func.(xs)
    display(fxs)

    sum(ws .* fxs)
end

export gauss_error
function gauss_error(exact, coefs, h;
    gauss_n = 5,
    x_begin=0, x_end=1,
    ux_begin=0, ux_end=0,
)
    N = length(coefs)

    phi1 = xi -> (1 - xi) / 2
    phi2 = xi -> (1 + xi) / 2

    x2xi = (xi, i) -> (h * ((1 + xi) / 2 + (i - 1))) * (x_end - x_begin) + x_begin

    ws, ps = gauss_quadrature_table[gauss_n]

    magic_coefs = cat(ux_begin, coefs, ux_end, dims=1)

    acc = 0.0

    for i = 1:(N+1)
        for j in 1:gauss_n
            diff = (
                exact(x2xi(ps[j], i))
                - magic_coefs[i] * phi1(ps[j])
                - magic_coefs[i+1] * phi2(ps[j])
            )
            acc += ws[j] * diff * diff
        end
    end
    acc *= h/2

    sqrt(acc)
end

end # module Common
