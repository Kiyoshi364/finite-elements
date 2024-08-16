module Hutils

export
    id,
    broadcast

id(x) = x
broadcast(f, xs) = f.(xs)

end # module Hutils
