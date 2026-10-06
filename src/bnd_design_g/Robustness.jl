using LinearAlgebra
using Statistics

"Modal-residue upper bound for ||C(sI-A)^(-1)B|| with C=I, B=scale*I."
function modal_residue_bound(A;scale=1.0)
    F=eigen(A); R=F.vectors
    try
        L=inv(R)'
        value=sum(norm(R[:,j])*norm(L[:,j])*scale for j in eachindex(F.values))
        return (;value,status="EVALUATED",condition=cond(R),eigenvalues=F.values)
    catch err
        return (;value=Inf,status="DEFECTIVE_OR_ILL_CONDITIONED",condition=Inf,
            eigenvalues=F.values,error=sprint(showerror,err))
    end
end

"Deterministically refine the largest singular value of a shifted resolvent."
function resolvent_peak(A,sigma,scale)
    λ=eigvals(A)
    alpha=maximum(real.(λ))
    if alpha >= -sigma-1e-11
        return (;status="NOMINAL_NOT_LEFT_OF_BOUNDARY",peak=Inf,omega_peak=NaN,
            beta_star=0.0,frequency=Float64[],sigma_min=Float64[],refinements=0)
    end
    As=A+sigma*I
    wmax=max(1.0,4maximum(abs.(imag.(λ)),init=0.0),20maximum(abs.(λ),init=0.0))
    grid=unique(vcat(0.0,10.0.^range(-4,log10(wmax),length=180),
        collect(range(0.0,wmax,length=240))))
    f(w)=scale*opnorm((im*w*I-As)\I,2)
    vals=f.(grid)
    peakix=argmax(vals); peak=vals[peakix]; wp=grid[peakix]; refined=0
    # Refine each sampled local maximum with deterministic golden-section steps.
    for i in 2:length(grid)-1
        if vals[i]>=vals[i-1] && vals[i]>=vals[i+1]
            a=grid[i-1]; b=grid[i+1]
            c=b-(b-a)*0.6180339887498949; d=a+(b-a)*0.6180339887498949
            fc=f(c); fd=f(d)
            for _ in 1:70
                if fc>fd
                    b=d; d=c; fd=fc; c=b-(b-a)*0.6180339887498949; fc=f(c)
                else
                    a=c; c=d; fc=fd; d=a+(b-a)*0.6180339887498949; fd=f(d)
                end
            end
            w=(a+b)/2; v=f(w); refined+=1
            if v>peak; peak=v; wp=w; end
        end
    end
    σmin=[1/max(v,eps()) for v in vals]
    (;status="EVALUATED",peak,omega_peak=wp,beta_star=1/peak,
      frequency=grid,sigma_min=σmin,refinements=refined)
end

"Full-complex-block small-gain certificate for the declared additive state-matrix ball."
function robustness_certificate(A,beta,sigma;scale=1.0)
    λ=eigvals(A)
    rp=resolvent_peak(A,sigma,scale)
    modal=modal_residue_bound(A;scale=scale)
    Rdelta=modal.value
    mrob=beta*Rdelta
    small_gain=isfinite(rp.peak) && beta*rp.peak<1
    (;status=rp.status,beta,scale,sigma_required=sigma,
        nominal_abscissa=maximum(real.(λ)),resolvent_peak=rp.peak,
        omega_peak=rp.omega_peak,beta_star=rp.beta_star,
        small_gain_pass=small_gain,small_gain_margin=1-beta*rp.peak,
        modal_residue_bound=Rdelta,modal_residue_status=modal.status,
        eigenvector_condition=modal.condition,m_rob=mrob,
        shifted_boundary_abscissa=maximum(real.(λ))+sigma,
        frequency=rp.frequency,sigma_min=rp.sigma_min,
        frequency_refinements=rp.refinements,
        certificate_model="A+scale*Delta, ||Delta||_2<=beta; complex full block")
end
