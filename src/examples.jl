module Examples

export Example, example
export bacarmo_example, bacarmo_example_gamma

struct Example
    f
    alpha :: Number
    beta :: Number
    gamma :: Number
    ux_begin :: Number
    ux_end :: Number
    x_begin :: Number
    x_end :: Number

    Example(f;
        alpha    :: Number = 1,
        beta     :: Number = 1,
        gamma    :: Number = 0,
        ux_begin :: Number = 0,
        ux_end   :: Number = 0,
        x_begin  :: Number = 0,
        x_end    :: Number = 1
    ) = new(
        f,
        alpha,
        beta,
        gamma,
        ux_begin,
        ux_end,
        x_begin,
        x_end,
    )

    Example(ex :: Example;
        f                                 = nothing,
        alpha    :: Union{Nothing,Number} = nothing,
        beta     :: Union{Nothing,Number} = nothing,
        gamma    :: Union{Nothing,Number} = nothing,
        ux_begin :: Union{Nothing,Number} = nothing,
        ux_end   :: Union{Nothing,Number} = nothing,
        x_begin  :: Union{Nothing,Number} = nothing,
        x_end    :: Union{Nothing,Number} = nothing
    ) = Example(
        (f     === nothing ? ex.f     : f    ),
        alpha   =(alpha    === nothing ? ex.alpha    : alpha   ),
        beta    =(beta     === nothing ? ex.beta     : beta    ),
        gamma   =(gamma    === nothing ? ex.gamma    : gamma   ),
        ux_begin=(ux_begin === nothing ? ex.ux_begin : ux_begin),
        ux_end  =(ux_end   === nothing ? ex.ux_end   : ux_end  ),
        x_begin =(x_begin  === nothing ? ex.x_begin  : x_begin ),
        x_end   =(x_end    === nothing ? ex.x_end    : x_end   ),
    )

    Base.broadcastable(x :: Example) = Ref(x)
end

struct Function
    exact
    deriv_2
end

function function_index(f_index :: UInt8, alpha :: Number, beta :: Number) :: Function
    i = f_index
    exact, deriv_2 = (
      i == 0 ? (x -> x + ((exp(-x) - exp(x)) / (exp(1) - exp(-1))), (x -> (exp(-x) - exp(x)) / (exp(1) - exp(-1)))) :
      i == 1 ? (x -> (-4) * x * (x - 1), (x -> -8)) :
      i == 2 ? ((x -> x * (x-1)), (x -> 2)) :
      i == 3 ? (x -> sin(pi * x), x -> - pi * pi * sin(pi * x)) :
      error("function_index out of bounds")
    )
    Function(exact, deriv_2)
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

function example(f_index :: UInt8, var_index :: UInt8 = 0) :: Tuple{Function, Example}
    ex = example_index(var_index)
    func = function_index(f_index, ex.alpha, ex.beta)
    func, Example(ex,
        ux_begin=func.exact(ex.x_begin),
        ux_end=func.exact(ex.x_end),
    )
end

# Examples from
# https://github.com/bacarmo/Problema-estacionario-unidimensional/blob/main/Eliptica_1D.ipynb
# in a different order
function bacarmo_example(index :: UInt8) :: Tuple{Any, Example}
    mk_f(alp, bet, exact, deriv_2) = x -> ((- alp) * deriv_2(x)) + (bet * exact(x))
    mk_ex(alp, bet, exact, deriv_2) = Example(mk_f(alp, bet, exact, deriv_2), alpha=alp, beta=bet)

    index == 0 ? begin
        alpha = 1.0
        beta = 1.0
        exact = x -> x + ((exp(-x) - exp(x)) / (exp(1) - exp(-1)))
        deriv_2 = x -> (exp(-x) - exp(x)) / (exp(1) - exp(-1))
        (exact, mk_ex(alpha, beta, exact, deriv_2))
    end : index == 1 ? begin
        alpha = 1.0
        beta = 0.0
        exact = x -> -4.0 * x * (x - 1.0)
        deriv_2 = x -> -8.0
        (exact, mk_ex(alpha, beta, exact, deriv_2))
    end : index == 2 ? begin
        alpha = 1.0
        beta = 1.0
        exact = x -> x * (x - 1.0)
        deriv_2 = x -> 2.0
        (exact, mk_ex(alpha, beta, exact, deriv_2))
    end : index == 3 ? begin
        alpha = 1.0
        beta = 1.0
        exact = x -> sin(pi * x)
        deriv_2 = x -> - pi * pi * sin(pi * x)
        (exact, mk_ex(alpha, beta, exact, deriv_2))
    end : error("function_index out of bounds")
end

# Examples from
# https://github.com/bacarmo/Problema-estacionario-unidimensional/blob/main/Eliptica_1D_caso2.ipynb
# in a different order
function bacarmo_example_gamma(index :: UInt8) :: Tuple{Any, Example}
    mk_f(alp, bet, gamm, exact, deriv_1, deriv_2) = x -> ((- alp) * deriv_2(x)) + (bet * exact(x)) + (gamm * deriv_1(x))
    mk_ex(alp, bet, gamm, exact, deriv_1, deriv_2) = Example(mk_f(alp, bet, gamm, exact, deriv_1, deriv_2), alpha=alp, beta=bet, gamma=gamm)

    index == 0 ? begin
        alpha = 1.0
        beta = 1.0
        gamma = 1.0
        exact = x -> x + ((exp(-x) - exp(x)) / (exp(1) - exp(-1)))
        deriv_1 = x -> 1 + ((- exp(-x)) - exp(x)) / (exp(1) - exp(-1))
        deriv_2 = x -> (exp(-x) - exp(x)) / (exp(1) - exp(-1))
        (exact, mk_ex(alpha, beta, gamma, exact, deriv_1, deriv_2))
    end : index == 1 ? begin
        alpha = 1.0
        beta = 0.0
        gamma = 1.0
        exact = x -> -4.0 * x * (x - 1.0)
        deriv_1 = x -> (-8.0 * x) + 4.0
        deriv_2 = x -> -8.0
        (exact, mk_ex(alpha, beta, gamma, exact, deriv_1, deriv_2))
    end : index == 2 ? begin
        alpha = 1.0
        beta = 1.0
        gamma = 2.0
        exact = x -> x * (x - 1.0)
        deriv_1 = x -> (2.0 * x) - 1.0
        deriv_2 = x -> 2.0
        (exact, mk_ex(alpha, beta, gamma, exact, deriv_1, deriv_2))
    end : index == 3 ? begin
        alpha = 1.0
        beta = 1.0
        gamma = 1.0
        exact = x -> sin(pi * x)
        deriv_1 = x -> pi * cos(pi * x)
        deriv_2 = x -> - pi * pi * sin(pi * x)
        (exact, mk_ex(alpha, beta, gamma, exact, deriv_1, deriv_2))
    end : error("function_index out of bounds")
end

end # module Examples

module TimeExamples

export Example, example
export bacarmo_example, bacarmo_example_gamma

struct Example
    f
    u0
    T :: Number
    alpha :: Number
    beta :: Number
    gamma :: Number

    Example(f, u0;
        T        :: Number = 1,
        alpha    :: Number = 1,
        beta     :: Number = 1,
        gamma    :: Number = 0,
    ) = new(
        f,
        u0,
        T,
        alpha,
        beta,
        gamma,
    )

    Example(ex :: Example;
        f                                 = nothing,
        u0                                = nothing,
        T        :: Union{Nothing,Number} = nothing,
        alpha    :: Union{Nothing,Number} = nothing,
        beta     :: Union{Nothing,Number} = nothing,
        gamma    :: Union{Nothing,Number} = nothing,
    ) = Example(
        (f     === nothing ? ex.f     : f    ),
        (u0    === nothing ? ex.u0    : u0   ),
        T       =(T        === nothing ? ex.T        : T       ),
        alpha   =(alpha    === nothing ? ex.alpha    : alpha   ),
        beta    =(beta     === nothing ? ex.beta     : beta    ),
        gamma   =(gamma    === nothing ? ex.gamma    : gamma   ),
    )

    Base.broadcastable(x :: Example) = Ref(x)
end

# Examples from
# https://github.com/bacarmo/Problema-estacionario-unidimensional/blob/main/Parabolica_1D_caso1.ipynb
# in a different order
function bacarmo_example(index :: UInt8) :: Tuple{Any, Example}
    mk_f(alp, bet, gamm, exact, dtime, deriv_1, deriv_2) = (x, t) -> ((- alp) * deriv_2(x, t)) + (bet * exact(x, t)) + (gamm * deriv_1(x, t) + dtime(x, t))
    mk_ex(_T, alp, bet, gamm, exact, dtime, deriv_1, deriv_2) = begin
        f = mk_f(alp, bet, 0.0, exact, dtime, deriv_1, deriv_2)
        u0 = x -> f(x, 0.0)
        Example(f, u0, alpha=alp, beta=bet, gamma=gamm, T=_T)
    end

    index == 0 ? begin
        T = 1.0
        alpha = 1.0
        beta = 1.0
        gamma = 0.0
        lamb = alpha * pi * pi + beta
        exact = (x, t) -> sin(pi * x) * (exp(lamb * t) / (pi * pi))
        dtime = (x, t) -> sin(pi * x) * lamb * (exp(lamb * t) / (pi * pi))
        deriv_1 = (x, t) -> pi * cos(pi * x) * (exp(lamb * t) / (pi * pi))
        deriv_2 = (x, t) -> - pi * pi * sin(pi * x) * (exp(lamb * t) / (pi * pi))
        (exact, mk_ex(alpha, beta, gamma, T, exact, dtime, deriv_1, deriv_2))
    end : index == 1 ? begin
        T = 1.0
        alpha = 1.0
        beta = 1.0
        gamma = 0.0
        lamb = 1.0
        exact = (x, t) -> sin(pi * x) * (exp(lamb * t) / (pi * pi))
        dtime = (x, t) -> sin(pi * x) * lamb * (exp(lamb * t) / (pi * pi))
        deriv_1 = (x, t) -> pi * cos(pi * x) * (exp(lamb * t) / (pi * pi))
        deriv_2 = (x, t) -> - pi * pi * sin(pi * x) * (exp(lamb * t) / (pi * pi))
        (exact, mk_ex(alpha, beta, gamma, T, exact, dtime, deriv_1, deriv_2))
    end : error("function_index out of bounds")
end

end # module TimeExamples
