module GraphReduction

using LinearAlgebra

export critical_graph_subspace, project_operator

function critical_graph_subspace(Phi,M,qvectors;threshold=0.99)
    isempty(qvectors) && throw(ArgumentError("at least one critical q vector is required"))
    energies=zeros(size(Phi,2))
    for q in qvectors
        z=Phi'*M*q; energies .+= abs2.(z)
    end
    energies ./= max(sum(energies),eps(Float64))
    order=sortperm(energies;rev=true)
    cumulative=cumsum(energies[order]); n=findfirst(>=(threshold),cumulative)
    selected=order[1:(n===nothing ? length(order) : n)]
    return (indices=selected,energy_fraction=sum(energies[selected]),mode_energies=energies)
end

project_operator(Phi_c,T) = Phi_c'*T*Phi_c

end
