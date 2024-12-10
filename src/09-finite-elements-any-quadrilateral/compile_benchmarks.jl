using BenchmarkTools: load
using BenchmarkTools: BenchmarkGroup
import BenchmarkTools as BT

using Plots: plot, plot!, savefig

const names = [
    :baseline, :ref, :iter, :iter_ref,
]
const markershapes = [
    :circle,
    :square,
    :diamond,
    :utriangle,
    :star5,
]
const name_prettyname = (;
    baseline = "Baseline",
    ref = "Ref",
    iter = "Iter",
    iter_ref = "Iter and Ref",
)

function make_graph(
    group :: BenchmarkGroup;
    label :: String = "label",
)
    local p = plot(
        xscale=:log2,
        xlabel="N1 (log2)",
        # yscale=:log10,
        ylabel="$(label)",
        legend=:topleft,
    )
    for (i, name) in enumerate(names)
        local subgroup = group[String(name)]
        local pretty_name = name_prettyname[name]
        local keys_str = [ k for k in keys(subgroup) ]
        local keys_int = parse.(Int64, keys_str)
        local perm = sortperm(keys_int)
        keys_str = keys_str[perm]
        keys_int = keys_int[perm]
        local values_ = getindex.((subgroup,), keys_str)
        plot!(p, keys_int, values_,
            label=pretty_name,
            markershape=markershapes[i],
        )
    end
    p
end

const bench = load("results.json")[1]

const file_label_func = [
    ("min_time", "minimum time (ns)", (x -> BT.time(BT.minimum(x)))),
    ("min_alloc", "allocs", (x -> BT.allocs(BT.minimum(x)))),
    ("min_memory", "memory (bytes)", (x -> BT.memory(BT.minimum(x)))),
]
const ext = "pdf"

for (file, label, f) in file_label_func
    local p = make_graph(
        f(bench),
        label=label
    )
    savefig(p, "$(file).$(ext)")
end
