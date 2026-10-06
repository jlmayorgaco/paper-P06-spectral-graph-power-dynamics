module ReplacementClosure

using LinearAlgebra

export port_resolvent, closure_matrix, determinant_lemma_residual

function port_resolvent(T0s,U,V)
    return V'*(T0s\U)
end

closure_matrix(Delta,G) = I+Delta*G

function determinant_lemma_residual(T0s,U,Delta,V)
    T=T0s+U*Delta*V'
    Tbase=T0s
    C=closure_matrix(Delta,port_resolvent(T0s,U,V))
    lt,pt=logabsdet(T); lb,pb=logabsdet(Tbase); lc,pc=logabsdet(C)
    log_res=abs(lt-lb-lc)
    phase_res=abs(pt-pb*pc)
    return (logabsdet_residual=log_res,phase_residual=phase_res,
            closure_dimension=size(C,1),T_sigma_min=minimum(svdvals(T)),
            C_sigma_min=minimum(svdvals(C)))
end

end
