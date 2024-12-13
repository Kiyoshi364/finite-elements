using BenchmarkTools: load
using BenchmarkTools: BenchmarkGroup
import BenchmarkTools as BT

using Plots: plot, plot!, savefig

const small_names = [
    :base, :base_precalc, :ref, :ref_precalc, :bacarmo
]
const smallname_prettyname = (;
    base = "Baseline",
    base_precalc = "Baseline with precalculation",
    ref = "Ref",
    ref_precalc = "Ref with precalculation",
    bacarmo = "Bruno Carmo",
)
const small_tests_names = [
    (:build_small_vec, "smallvec"),
    (:build_small_mat, "smallmat"),
]
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

function make_small_csv(
    group :: BenchmarkGroup,
    io :: IO,
)
    local colunms_funcs = [
        ("min_time (ns)", (x -> BT.time(BT.minimum(x)))),
        ("allocs", (x -> BT.allocs(BT.minimum(x)))),
        ("memory (bytes)", (x -> BT.memory(BT.minimum(x)))),
    ]
    write(io, "Name")
    for (col_name, func) in colunms_funcs
        write(io, ",")
        write(io, col_name)
    end
    write(io, "\n")
    for small_name in small_names
        local trial = group[String(small_name)]
        write(io, smallname_prettyname[small_name])
        for (col_name, func) in colunms_funcs
            write(io, ",")
            write(io, string(func(trial)))
        end
        write(io, "\n")
    end
end

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

for (small_test, small_test_name) in small_tests_names
    local small_test_bench = bench[small_test]
    if 0 < length(small_test_bench)
        local file = open(
            "$(out_dir)$(small_test_name).csv",
            lock = false,
            create = true,
            write = true,
        )
        make_small_csv(small_test_bench, file)
        close(file)
    end
end

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
