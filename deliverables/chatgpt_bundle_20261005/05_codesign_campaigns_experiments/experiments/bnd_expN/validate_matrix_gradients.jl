using CSV, DataFrames, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
isfile(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json")) ||
    error("post-freeze derivative addendum requires frozen model")
ctx=design_context(ROOT)
rows=NamedTuple[]
for name in ("random_1","random_2","random_3")
    t=CSV.read(joinpath(ROOT,"reports","experiment_N","validation_cases",
        name,"case_definition.csv"),DataFrame)
    rho=Float64.(t.rho);kp=Float64.(t.Kp);ki=Float64.(t.Ki)
    sp=spectrum(ctx,rho,kp,ki)
    d=simple_mode_sensitivities(ctx,rho,kp,ki)
    Q=sp.quotient;Aq=transpose(Q)*sp.model.Ared*Q
    F=eigen(Aq);j=argmin(abs.(F.values.-d.lambda))
    right=F.vectors[:,j]
    FL=eigen(adjoint(Aq));left=FL.vectors[:,argmin(abs.(FL.values.-conj(d.lambda)))]
    den=dot(left,right)
    for (kind,values,exact) in (("rho",rho,d.rho),("Kp",kp,d.Kp),("Ki",ki,d.Ki))
        for i in 1:10
            kind=="rho" && i!=9 && !(0<rho[i]<1) && continue
            h=kind=="rho" ? 1e-4 : 1e-3*max(1.0,abs(values[i]))
            plus=copy(values);minus=copy(values)
            plus[i]+=h;minus[i]-=h
            mp=kind=="rho" ? descriptor(ctx,plus,kp,ki) :
                kind=="Kp" ? descriptor(ctx,rho,plus,ki) :
                             descriptor(ctx,rho,kp,plus)
            mm=kind=="rho" ? descriptor(ctx,minus,kp,ki) :
                kind=="Kp" ? descriptor(ctx,rho,minus,ki) :
                             descriptor(ctx,rho,kp,minus)
            dA=(mp.Ared-mm.Ared)/(2h)
            fd=dot(left,transpose(Q)*dA*Q*right)/den
            err=abs(fd-exact[i])
            rel=err/max(abs(exact[i]),1e-8)
            push!(rows,(;case=name,bus=i+29,parameter=kind,step=h,
                exact_real=real(exact[i]),exact_imag=imag(exact[i]),
                matrix_FD_real=real(fd),matrix_FD_imag=imag(fd),
                absolute_error=err,relative_error=rel,
                status=(rel<1e-5 || err<1e-7) ? "PASS" : "FAIL"))
        end
    end
    println("MATRIX_GRADIENT_CASE ",name);flush(stdout)
end
CSV.write(joinpath(ROOT,"reports","experiment_N",
    "TABLE_N08_matrix_derivative_addendum.csv"),DataFrame(rows))
println("MATRIX_GRADIENT_DONE rows=",length(rows)," failures=",
    count(r->r.status!="PASS",rows)," max_abs=",
    maximum(r.absolute_error for r in rows)," max_rel=",
    maximum(r.relative_error for r in rows));flush(stdout)
