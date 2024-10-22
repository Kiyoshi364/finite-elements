module Examples2d

export Example
export example

struct Example
    f
    alpha :: Number
    beta :: Number

    Example(f;
        alpha    :: Number = 1,
        beta     :: Number = 1,
    ) = new(
        f,
        alpha,
        beta,
    )

    Example(ex :: Example;
        f                                 = nothing,
        alpha    :: Union{Nothing,Number} = nothing,
        beta     :: Union{Nothing,Number} = nothing,
    ) = Example(
        (f         === nothing ? ex.f          : f    ),
        alpha   =(alpha    === nothing ? ex.alpha    : alpha   ),
        beta    =(beta     === nothing ? ex.beta     : beta    ),
    )

    Base.broadcastable(x :: Example) = Ref(x)
end

function example(index :: UInt8) :: Tuple{Any, Example}
    mk_f(alp, bet, exact, deriv_2) =
        (x) -> (
            ((- alp) * foldl(+, map(i -> deriv_2[i](x), 1:length(x))))
            + (bet * exact(x))
        )
    mk_ex(alp, bet, exact, deriv_2) = begin
        f = mk_f(alp, bet, exact, deriv_2)
        Example(f, alpha=alp, beta=bet)
    end

    index == 0 ? begin
        alpha = 1.0
        beta = 1.0
        exact = x -> sin(pi * x[1]) * sin(pi * x[2])
        deriv_1 = [
            x -> pi * cos(pi * x[1]) * sin(pi * x[2]),
            x -> pi * sin(pi * x[1]) * cos(pi * x[2]),
        ]
        deriv_2 = [
            x -> (-1 * pi * pi) * sin(pi * x[1]) * sin(pi * x[2]),
            x -> (-1 * pi * pi) * sin(pi * x[1]) * sin(pi * x[2]),
        ]
        (exact, mk_ex(alpha, beta, exact, deriv_2))
    # end : index == 1 ? begin
    end : error("function_index out of bounds")
end

end # module Examples2d
