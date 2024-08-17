module FiniteDifferences

export finite_differences
export Example, example

using LinearAlgebra

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

function build_mat(alpha, beta, hsqr, dim)
    #    (,~ $ ] }.@,@# 1 2 1 ,:@, 0 $~ 0>.-&2)
    @assert dim >= 2
    @assert alpha > 0
    a = - alpha
    b = (2 * alpha) + (beta * hsqr)
    zeros = fill(0, (dim-2,))
    prefix = [a, b, a]

    template = cat(prefix, zeros, dims=1)
    templates = fill(template, (dim,))
    flat = cat(templates..., dims=1)
    un_headed = flat[2:end-(dim-1)]
    ret = reshape(un_headed, (dim, dim))
    ret
end

function build_vec(f, hsqr, dim;
    alpha=1,
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert dim >= 2

    offset = cat(
        [ux_begin * alpha],
        fill(0, (dim-2,)),
        [ux_end * alpha],
    dims=1)

    inter = (1:dim) ./ (dim+1)
    xs = (x_end - x_begin) * inter .+ x_begin
    fxs = f.(xs)
    fs = fxs * hsqr
    sum = fs + offset

    ret = sum
    ret
end

function finite_differences(ex :: Example, h, N)
    finite_differences(ex.f, ex.alpha, ex.beta, h, N;
        ux_begin=ex.ux_begin, ux_end=ex.ux_end,
        x_begin=ex.x_begin, x_end=ex.x_end,
    )
end

function finite_differences(f, alpha, beta, h, N;
    ux_begin=0, ux_end=0,
    x_begin=0, x_end=1,
)
    @assert h == (x_end - x_begin) / (N + 1)
    hsqr = h * h
    A = build_mat(alpha, beta, hsqr, N)
    b = build_vec(f, hsqr, N,
        alpha=alpha,
        ux_begin=ux_begin, ux_end=ux_end,
        x_begin=x_begin, x_end=x_end,
    )
    uh = A \ b
    uh
end

function function_index(f_index :: UInt8, alpha :: Number, beta :: Number)
    i = f_index
    exact, deriv_2 = (
      i == 0 ? ((x -> x * (x-1)), (x -> 2)) :
      i == 1 ? (x -> sin(pi * x), x -> - pi * pi * sin(pi * x)) :
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

end # module FiniteDifferences
