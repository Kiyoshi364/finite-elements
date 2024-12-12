using BenchmarkTools: load
using BenchmarkTools: BenchmarkGroup
import BenchmarkTools as BT

using Plots: plot, plot!, savefig

const names = [
    :baseline, :ref, :iter, :iter_ref, :bacarmo
]
const tests_names = [
    (:build_vec, "vec"),
    (:build_mat, "mat"),
    (:build_both, "both"),
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
    bacarmo = "Bruno Carmo",
)

function make_group_graph(
    group :: BenchmarkGroup;
    label :: String = "label",
)
    local scale = :log10
    local p = plot(
        xscale=scale,
        xlabel="N ($(String(scale)))",
        yscale=scale,
        ylabel="$(label) ($(String(scale)))",
        legend=:topleft,
    )
    for (i, name) in enumerate(names)
        local subgroup = group[String(name)]
        if 0 < length(subgroup)
            local pretty_name = name_prettyname[name]
            local keys_str = [ k for k in keys(subgroup) ]
            local keys_int = parse.(Int64, keys_str)
            local perm = sortperm(keys_int)
            keys_str = keys_str[perm]
            keys_int = keys_int[perm]
            local values_ = getindex.((subgroup,), keys_str)
            plot!(p, keys_int .* keys_int, values_,
                label=pretty_name,
                markershape=markershapes[i],
            )
        end
    end
    p
end

function make_single_graph(
    name :: String,
    group :: BenchmarkGroup,
)
    local label_func = [
        ("minimum time (ns)", (x -> BT.time(BT.minimum(x)))),
        ("allocs", (x -> BT.allocs(BT.minimum(x)))),
        ("memory (bytes)", (x -> BT.memory(BT.minimum(x)))),
    ]
    local scale = :log10
    local p = plot(
        xscale=scale,
        xlabel="N elements ($(String(scale)))",
        yscale=scale,
        ylabel="($(String(scale)))",
        legend=:topleft,
    )
    local keys_str = [ k for k in keys(group) ]
    local keys_int = parse.(Int64, keys_str)
    local perm = sortperm(keys_int)
    keys_str = keys_str[perm]
    keys_int = keys_int[perm]
    local trials = getindex.((group,), keys_str)
    for (i, (label, func)) in enumerate(label_func)
        local values_ = func.(trials)
        plot!(p, keys_int .* keys_int, values_,
            label=label,
            markershape=markershapes[i],
        )
    end
    p
end

const bench = load("results.json")[1]
const out_dir = "./out/"
const ext = "pdf"

const groups_file_label_func = [
    ("min_time", "minimum time (ns)", (x -> BT.time(BT.minimum(x)))),
    ("min_alloc", "allocs", (x -> BT.allocs(BT.minimum(x)))),
    ("min_memory", "memory (bytes)", (x -> BT.memory(BT.minimum(x)))),
]

for (test, test_name) in tests_names
    local test_bench = bench[test]
    for (file, label, f) in groups_file_label_func
        local p = make_group_graph(
            f(test_bench),
            label=label
        )
        savefig(p, "$(out_dir)$(test_name)-$(file).$(ext)")
    end

    for (name, single_bench) in test_bench
        if 0 < length(single_bench)
            local p = make_single_graph(
                name,
                single_bench,
            )
            savefig(p, "$(out_dir)$(test_name)-$(name).$(ext)")
        end
    end
end
