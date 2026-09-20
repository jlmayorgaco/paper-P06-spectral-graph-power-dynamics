"""Locate the campaign root in the source tree or a freshly extracted bundle."""

function campaign_root(; explicit=nothing)
    if explicit !== nothing
        root = abspath(String(explicit))
        isdir(joinpath(root, "raw")) && isdir(joinpath(root, "reports")) ||
            error("not an IAS2026 campaign root: $root")
        return normpath(root)
    end
    start = abspath(@__DIR__)
    current = start
    while true
        if isdir(joinpath(current, "raw")) && isdir(joinpath(current, "reports")) &&
           (isdir(joinpath(current, "src")) || isdir(joinpath(current, "code")))
            return normpath(current)
        end
        parent = dirname(current)
        parent == current && break
        current = parent
    end
    error("could not locate campaign root from Julia script path")
end

function campaign_root_from_args()
    idx = findfirst(==("--campaign-root"), ARGS)
    idx === nothing && return campaign_root()
    idx < length(ARGS) || error("--campaign-root requires a path")
    return campaign_root(explicit=ARGS[idx + 1])
end
