using LinearAlgebra,CSV,DataFrames,Statistics
include(joinpath(@__DIR__,"DDE_NEV.jl"))
const N=DDENEV;const D=N.D;const R=D.R;const BASE=joinpath(@__DIR__,"baseline_reproduction")
BLAS.set_num_threads(1);ctx=R.N.design_context(R.ROOT)
rho=fill(.875,10);kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
frozen=CSV.read(joinpath(@__DIR__,"FROZEN_DELAY_PATTERNS.csv"),DataFrame)
cases=[("uniform_10ms",fill(10.0,10)),
       ("min_commutator",parse.(Float64,split(only(frozen.tau_ms[frozen.delay_pattern_id.=="min_commutator"]),";"))),
       ("max_commutator",parse.(Float64,split(only(frozen.tau_ms[frozen.delay_pattern_id.=="max_commutator"]),";")))]
buses=[30,32,33];rows=NamedTuple[]
for (case,tms) in cases
    tau=tms./1000;L=D.linearization(ctx,rho,kp,ki)
    roots=N.track_roots(L,tau;real_cut=-.25,continuation_steps=4,max_roots=2,maxiter=12)
    root=only(filter(r->imag(r.s)>0,roots));M=D.delta_matrix(L,root.s,tau);S=svd(M)
    w=S.U[:,end];v=root.v;denom=dot(w,D.delta_s(L,root.s,tau)*v)
    for bus in buses
        i=bus-29
        derivatives=Dict(
            "tau"=>D.delta_tau(L,root.s,tau,i),
            "Kp"=>(-exp(-root.s*tau[i])).*(L.Bp[:,i]*L.C[:,i]'),
            "Ki"=>(-exp(-root.s*tau[i])).*(L.Bi[:,i]*L.C[:,i]'))
        params=[("tau",tau[i],2e-4),("Kp",kp[i],1e-2*kp[i]),("Ki",ki[i],1e-2*ki[i])]
        for (parameter,pvalue,step) in params
            analytic=-dot(w,derivatives[parameter]*v)/denom
            if parameter=="tau"
                tp=copy(tau);tm=copy(tau);tp[i]+=step;tm[i]-=step
                rp=N.newton_root(L,root.s,tp,root.v;tol=1e-14,maxiter=20);rm=N.newton_root(L,root.s,tm,root.v;tol=1e-14,maxiter=20)
            else
                kpp=copy(kp);kpm=copy(kp);kip=copy(ki);kim=copy(ki)
                if parameter=="Kp";kpp[i]+=step;kpm[i]-=step;else;kip[i]+=step;kim[i]-=step;end
                Lp=D.linearization(ctx,rho,kpp,kip);Lm=D.linearization(ctx,rho,kpm,kim)
                rp=N.newton_root(Lp,root.s,tau,root.v;tol=1e-14,maxiter=20);rm=N.newton_root(Lm,root.s,tau,root.v;tol=1e-14,maxiter=20)
            end
            converged=rp.converged&&rm.converged
            fd=converged ? (rp.s-rm.s)/(2step) : complex(NaN,NaN)
            abserr=converged ? abs(analytic-fd) : Inf
            relerr=converged ? abserr/max(abs(analytic),abs(fd),1e-10) : Inf
            push!(rows,(;delay_pattern_id=case,bus,parameter,value=pvalue,fd_step=step,
                tracked_root_real=real(root.s),tracked_root_imag=imag(root.s),
                analytic_real=real(analytic),analytic_imag=imag(analytic),fd_real=real(fd),fd_imag=imag(fd),
                absolute_error=abserr,relative_error=relerr,left_right_denominator=abs(denom),
                plus_minus_roots_converged=converged,root_nev_residual=root.residual,
                sensitivity_status=converged&&isfinite(relerr)&&relerr<=.01 ? "NUMERICALLY_VALIDATED" : "SUPPORTED_LOCAL_OR_FAILED"))
        end
    end
    println("D02_LOCAL_SENSITIVITY_DONE ",case);flush(stdout)
end
CSV.write(joinpath(BASE,"TABLE_D02_DERIVATIVE_VALIDATION.csv"),DataFrame(rows))
status=(;validated_count=count(r->r.sensitivity_status=="NUMERICALLY_VALIDATED",rows),
    attempted_count=length(rows),parameter_coverage="Kp, Ki, tau only; rho sensitivity not derived/validated",
    overall_status="BLOCKED_PARAMETER_SET_INCOMPLETE",root_scope="two locally tracked roots; no DDE spectral completeness")
open(joinpath(BASE,"TABLE_D02_STATUS.toml"),"w") do io
    for (k,v) in pairs(status);println(io,k," = ",v isa AbstractString ? "\"$v\"" : v);end
end
println(status)
