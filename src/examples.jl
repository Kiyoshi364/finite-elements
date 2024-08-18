module Examples

export Example, example

struct Example
    f
    alpha :: Number
    beta :: Number
    ux_begin :: Number
    ux_end :: Number
    x_begin :: Number
    x_end :: Number

    Example(f;
        alpha    :: Number = 1,
        beta     :: Number = 1,
        ux_begin :: Number = 0,
        ux_end   :: Number = 0,
        x_begin  :: Number = 0,
        x_end    :: Number = 1
    ) = new(
        f,
        alpha,
        beta,
        ux_begin,
        ux_end,
        x_begin,
        x_end,
    )

    Example(ex :: Example;
        f                                 = nothing,
        alpha    :: Union{Nothing,Number} = nothing,
        beta     :: Union{Nothing,Number} = nothing,
        ux_begin :: Union{Nothing,Number} = nothing,
        ux_end   :: Union{Nothing,Number} = nothing,
        x_begin  :: Union{Nothing,Number} = nothing,
        x_end    :: Union{Nothing,Number} = nothing
    ) = Example(
        (f     === nothing ? ex.f     : f    ),
        alpha   =(alpha    === nothing ? ex.alpha    : alpha   ),
        beta    =(beta     === nothing ? ex.beta     : beta    ),
        ux_begin=(ux_begin === nothing ? ex.ux_begin : ux_begin),
        ux_end  =(ux_end   === nothing ? ex.ux_end   : ux_end  ),
        x_begin =(x_begin  === nothing ? ex.x_begin  : x_begin ),
        x_end   =(x_end    === nothing ? ex.x_end    : x_end   ),
    )

    Base.broadcastable(x :: Example) = Ref(x)
end

function function_index(f_index :: UInt8, alpha :: Number, beta :: Number)
    i = f_index
    exact, deriv_2 = (
      i == 0 ? (x -> x + ((exp(-x) - exp(x)) / (exp(1) - exp(-1))), (x -> (exp(-x) - exp(x)) / (exp(1) - exp(-1)))) :
      i == 1 ? ((x -> x * (x-1)), (x -> 2)) :
      i == 2 ? (x -> sin(pi * x), x -> - pi * pi * sin(pi * x)) :
      error("function_index out of bounds")
    )
    f = x -> ((- alpha) * deriv_2(x)) + (beta * exact(x))
    exact, f
end

function example_index(var_index :: UInt8) :: Example
    bsize = 0

    alpha = begin
        alloc = 1
        i = (var_index >> bsize) & ((1 << alloc) - 1)
        bsize += alloc
        (
            i == 0 ? 1.0 :
            i == 1 ? 7.0 :
            error("alpha index out of bounds")
        )
    end

    beta = begin
        alloc = 1
        i = (var_index >> bsize) & ((1 << alloc) - 1)
        bsize += alloc
        (
            i == 0 ? 1.0 :
            i == 1 ? 0.0 :
            error("beta index out of bounds")
        )
    end

    x_begin = begin
        alloc = 2
        i = (var_index >> bsize) & ((1 << alloc) - 1)
        bsize += alloc
        (
            i == 0 ? 0.0 :
            i == 1 ? 2.0 :
            i == 2 ? -1.0 :
            i == 3 ? -2.0 :
            error("x_begin index out of bounds")
        )
    end

    x_end = begin
        alloc = 2
        i = (var_index >> bsize) & ((1 << alloc) - 1)
        bsize += alloc
        (
            i == 0 ? 1.0 :
            i == 1 ? 2.5 :
            i == 2 ? 0.5 :
            i == 3 ? -1.5 :
            error("x_end index out of bounds")
        )
    end

    Example(nothing,
        alpha=alpha,
        beta=beta,
        x_begin=x_begin,
        x_end=x_end,
    )
end

function example(f_index :: UInt8, var_index :: UInt8 = 0) :: Tuple{Any, Example}
    ex = example_index(var_index)
    exact, f = function_index(f_index, ex.alpha, ex.beta)
    exact, Example(ex, f=f,
        ux_begin=exact(ex.x_begin),
        ux_end=exact(ex.x_end),
    )
end

end # module Examples
