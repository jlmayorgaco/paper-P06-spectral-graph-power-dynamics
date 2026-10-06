module IdealDamping

using LinearAlgebra

export modal_critical_damping, modal_poles, isolated_decay_rate

function symmetric_power(A, power; positive_semidefinite=false, tol=1e-11)
    E=eigen(Hermitian((A+A')/2))
    scale=max(maximum(abs,E.values;init=0.0),eps(Float64))
    minimum(E.values)>=-tol*scale || throw(ArgumentError("matrix has negative eigenvalues"))
    if positive_semidefinite
        vals=max.(E.values,0.0).^power
    else
        minimum(E.values)>tol*scale || throw(ArgumentError("matrix must be positive definite"))
        vals=E.values.^power
    end
    return E.vectors*Diagonal(vals)*E.vectors'
end

"""Return modal critical damping for M>0, K>=0 without a global-optimum claim."""
function modal_critical_damping(M::AbstractMatrix,K::AbstractMatrix;tol=1e-11)
    size(M)==size(K) || throw(DimensionMismatch("M and K must have same size"))
    Mhalf=symmetric_power(M,0.5;tol)
    Minvhalf=symmetric_power(M,-0.5;tol)
    Ktilde=Minvhalf*K*Minvhalf
    Ktilde=Hermitian((Ktilde+Ktilde')/2)
    E=eigen(Ktilde)
    scale=max(maximum(abs,E.values;init=0.0),eps(Float64))
    minimum(E.values)>=-tol*scale || throw(ArgumentError("K must be positive semidefinite"))
    ν=max.(E.values,0.0)
    Khalf=E.vectors*Diagonal(sqrt.(ν))*E.vectors'
    D=2*Mhalf*Khalf*Mhalf
    return (D=Matrix(Hermitian((D+D')/2)),nu=ν,U=E.vectors,
        Mhalf=Mhalf,Minvhalf=Minvhalf,Ktilde=Matrix(Ktilde))
end

function modal_poles(nu,d)
    nu>=0 && d>=0 || throw(ArgumentError("nu and d must be nonnegative"))
    disc=complex(d^2-4nu)
    return ((-d+sqrt(disc))/2,(-d-sqrt(disc))/2)
end

"""Asymptotic decay rate of the slow root of s^2+d*s+nu=0."""
function isolated_decay_rate(nu,d)
    λ1,λ2=modal_poles(nu,d)
    return -maximum(real.((λ1,λ2)))
end

end
