module RealPoleAnalysis

using LinearAlgebra
using Statistics

export pole_metrics, gauge_score, select_real_poles, graph_projection,
       spearman_rank

pole_metrics(s::Number) = (frequency_hz=abs(imag(s))/(2pi),
    damping_ratio=abs(s) == 0 ? NaN : -real(s)/abs(s))

function gauge_score(q::AbstractVector, g::AbstractVector)
    return abs(dot(g, q))^2 / max(norm(g)^2*norm(q)^2, eps(Float64))
end

function graph_projection(Phi, M, q::AbstractVector)
    z = Phi' * M * q
    p = abs2.(z)
    total = sum(p)
    total > 0 || return (coefficients=z, fractions=fill(NaN,length(z)),
                        concentration=NaN, top_mode=0, top1_top2=NaN)
    p ./= total
    order = sortperm(p; rev=true)
    ratio = length(order) < 2 ? Inf : p[order[1]]/max(p[order[2]],eps(Float64))
    return (coefficients=z, fractions=p, concentration=maximum(p),
            top_mode=order[1], top1_top2=ratio)
end

"""Select non-conjugate actual poles by preregistered spectral roles."""
function select_real_poles(A::AbstractMatrix, state_names::Vector{String},
                           qidx::Vector{Int}, vidx::Vector{Int},
                           Phi::AbstractMatrix, M::AbstractMatrix;
                           low_hz=0.05, high_hz=5.0, keep_band=4)
    E = eigen(A)
    representatives = findall(i -> imag(E.values[i]) > 1e-7, eachindex(E.values))
    rows = NamedTuple[]
    n = size(Phi,1)
    g = ones(n)
    for i in representatives
        s = E.values[i]
        x = E.vectors[:,i]
        q = x[qidx]
        v = x[vidx]
        pm = pole_metrics(s)
        gp = graph_projection(Phi,M,q)
        retained = sum(abs2,x[vcat(qidx,vidx)]) / max(sum(abs2,x),eps(Float64))
        weights = abs2.(x)
        dom = sortperm(weights;rev=true)[1:min(3,length(weights))]
        names = [state_names[j] for j in dom]
        pllidx = findall(nm -> occursin("pll₊",nm) || occursin("pll+",nm), state_names)
        pll_share = isempty(pllidx) ? 0.0 : sum(weights[intersect(pllidx, eachindex(weights))]) / max(sum(weights),eps(Float64))
        push!(rows,(index=i,lambda=s,frequency_hz=pm.frequency_hz,
            damping_ratio=pm.damping_ratio,retained_participation=retained,
            gauge_score=gauge_score(q,g),q_norm=norm(q),projection=gp.fractions,
            dominant_graph_mode=gp.top_mode,graph_concentration=gp.concentration,
            top1_top2=gp.top1_top2,pll_state_share=pll_share,
            dominant_states=join(names,"; "),right_vector=x))
    end
    selected = Dict{Int,Tuple{String,NamedTuple}}()
    function add_selected(r, role)
        if haskey(selected,r.index)
            oldrole,oldrow=selected[r.index]
            roles=split(oldrole,"; ")
            role in roles || push!(roles,role)
            selected[r.index]=(join(roles,"; "),oldrow)
        else
            selected[r.index]=(role,r)
        end
    end
    non_gauge = filter(r -> r.gauge_score < 0.9 && abs(r.lambda)>1e-7, rows)
    if !isempty(non_gauge)
        sp = non_gauge[argmax(real.(getproperty.(non_gauge,:lambda)))]
        add_selected(sp,"spectral_abscissa_pair")
    end
    band = filter(r -> low_hz <= r.frequency_hz <= high_hz &&
                       r.gauge_score < 0.9 && r.retained_participation > 1e-6, rows)
    sort!(band;by=r -> (r.damping_ratio,-r.retained_participation))
    for r in band[1:min(keep_band,length(band))]
        add_selected(r,"lowest_damping_band_pair")
    end
    remaining = filter(r -> r.gauge_score < 0.9 && r.retained_participation > 1e-3 &&
                            !haskey(selected,r.index), rows)
    sort!(remaining;by=r -> -r.retained_participation)
    for r in remaining[1:min(2,length(remaining))]
        add_selected(r,"high_retained_participation")
    end
    pll = filter(r -> r.gauge_score < 0.9 && r.pll_state_share > 0, rows)
    if !isempty(pll)
        sort!(pll;by=r -> -r.pll_state_share)
        r=first(pll)
        add_selected(r,"PLL_state_dominant")
    end
    result = NamedTuple[]
    for (idx,(role,r)) in selected
        push!(result, merge(r,(role=role,)) )
    end
    sort!(result;by=r -> r.frequency_hz)
    return (all_positive=rows, selected=result, values=E.values, vectors=E.vectors)
end

function spearman_rank(x::AbstractVector, y::AbstractVector)
    length(x)==length(y) || throw(ArgumentError("rank vectors must match"))
    length(x)<2 && return NaN
    function ranks(v)
        p=sortperm(v)
        r=zeros(Float64,length(v))
        i=1
        while i<=length(v)
            j=i
            while j<length(v) && v[p[j+1]]==v[p[i]]; j+=1; end
            r[p[i:j]].=(i+j)/2
            i=j+1
        end
        return r
    end
    rx,ry=ranks(x),ranks(y)
    return cor(rx,ry)
end

end
