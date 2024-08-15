using LinearAlgebra

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
    return ret
end

function build_vec(f, alpha, a, b, hsqr, dim)
    @assert dim >= 2

    offset = cat(
        [f(a) * alpha],
        fill(0, (dim-2,)),
        [f(b) * alpha],
    dims=1)

    inter = (1:dim) ./ (dim+1)
    xs = (b - a) * inter .+ a
    fxs = f.(xs)
    fs = fxs * hsqr
    sum = fs .+ offset

    ret = sum
    return ret
end

bgn, nd = 0, 1
alpha = 1
beta = 2

if true
    exact(x) = x * (x-1)
    deriv_2(x) = -2
else
    exact(x) = sin(pi * x)
    deriv_2(x) = pi * pi * sin(pi * x)
end
func(x) = (deriv_2(x) * alpha) + (beta * exact(x))

last_err = nothing
for i = 3:10
    N = (1 << i) - 1
    h = (nd - bgn) / (N + 1)
    hsqr = h*h
    A = build_mat(alpha, beta, hsqr, N)
    # display(A)
    b = build_vec(func, alpha, bgn, nd, hsqr, N)
    # display(b)
    uh = A \ b
    # display(uh)
    u = build_vec(exact, 0, bgn, nd, 1, N)
    # display(u)
    # err = norm(uh - u)/norm(u)
    err = maximum(abs.(uh - u))
    if last_err == nothing
        nothing
    else
        print("div: ")
        println(last_err/err)
    end
    global last_err = err
    println(err)
end
