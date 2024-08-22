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

end # module Common
