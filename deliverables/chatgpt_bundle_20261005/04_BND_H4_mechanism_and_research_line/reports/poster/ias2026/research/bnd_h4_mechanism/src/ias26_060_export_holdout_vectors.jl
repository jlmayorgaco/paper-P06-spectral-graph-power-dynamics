using CSV
using DataFrames
using JLD2
using LinearAlgebra

function export_holdout(runroot::String)
    selection = CSV.read(joinpath(runroot, "tables", "CROSSCODE_HOLDOUT_SELECTION.csv"), DataFrame)
    selected = selection[selection.selection_status .!= "EMPTY_STRATUM_NO_REPLACEMENT", :]
    ids = unique(String.(selected.scenario_id))
    length(ids) == 50 || error("expected 50 unique holdout IDs, got $(length(ids))")
    outdir = joinpath(runroot, "raw", "crosscode")
    mkpath(outdir)
    outpath = joinpath(outdir, "JULIA_CRITICAL_EIGENVECTORS.csv")
    rows = NamedTuple[]
    for sid in ids
        path = joinpath(runroot, "raw", "julia", "scenarios", "$(sid).jld2")
        isfile(path) || error("missing immutable Julia spectrum for $sid")
        jldopen(path, "r") do f
            String(f["metadata/scenario_id"]) == sid || error("scenario ID mismatch in $path")
            prefixes = String.(f["metadata/case_prefixes"])
            for (treatment, idx) in (("ORIGINAL", 16), ("INTERVENTION", 17))
                prefix = prefixes[idx]
                startswith(prefix, "cases/$(lpad(string(idx), 2, '0'))_") || error("unexpected portfolio ordering: $prefix")
                values = ComplexF64.(f["$prefix/eigenvalues"])
                vectors = ComplexF64.(f["$prefix/eigenvectors"])
                mask = findall(abs.(values) .>= 1e-3)
                isempty(mask) && error("empty transverse spectrum for $sid $treatment")
                critical = mask[argmax(real.(values[mask]))]
                vector = vectors[:, critical]
                lam = values[critical]
                for j in eachindex(vector)
                    push!(rows, (scenario_id=sid, treatment=treatment, vector_index=j,
                        state_count=length(vector), lambda_real=real(lam), lambda_imag=imag(lam),
                        vector_real=real(vector[j]), vector_imag=imag(vector[j])))
                end
            end
        end
    end
    CSV.write(outpath, DataFrame(rows))
    println("IAS26_060_JULIA_HOLDOUT_VECTORS ids=$(length(ids)) rows=$(length(rows)) path=$outpath")
end

length(ARGS) == 1 || error("usage: julia --project=<env> ias26_060_export_holdout_vectors.jl <RUN_ROOT>")
export_holdout(abspath(ARGS[1]))
