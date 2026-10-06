module PowerDynamicsValidation

using SHA
using CSV
using DataFrames
using LinearAlgebra

export candidate_gate, run_validation

sha256_file(path)=bytes2hex(sha256(read(path)))

"Check that the analytic candidate is present, frozen, and matches its sidecar hash."
function candidate_gate(root)
    candidate=joinpath(root,"reports","experiment_D","ANALYTIC_CANDIDATE.json")
    sidecar=joinpath(root,"reports","experiment_D","ANALYTIC_CANDIDATE.sha256")
    isfile(candidate) || return (allowed=false,reason="candidate file is absent",candidate=candidate)
    isfile(sidecar) || return (allowed=false,reason="candidate SHA-256 freeze record is absent",candidate=candidate)
    raw=read(candidate,String)
    occursin(r"\"candidate_frozen\"\s*:\s*true",raw) ||
        return (allowed=false,reason="candidate_frozen is not true",candidate=candidate)
    occursin(r"\"status\"\s*:\s*\"FROZEN\"",raw) ||
        return (allowed=false,reason="candidate status is not FROZEN",candidate=candidate)
    expected=lowercase(strip(read(sidecar,String)))
    actual=sha256_file(candidate)
    expected==actual || return (allowed=false,reason="candidate hash differs from freeze record",candidate=candidate)
    return (allowed=true,reason="frozen candidate hash verified",candidate=candidate,sha256=actual,raw=raw)
end

function json_number(object::AbstractString,key::AbstractString)
    pattern=Regex("\\\""*key*"\\\"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)")
    m=match(pattern,object)
    m===nothing && error("frozen candidate is missing numeric field $key")
    return parse(Float64,m.captures[1])
end

function json_bool(object::AbstractString,key::AbstractString)
    m=match(Regex("\\\""*key*"\\\"\\s*:\\s*(true|false)"),object)
    m===nothing && error("frozen candidate is missing boolean field $key")
    return m.captures[1]=="true"
end

function candidate_rows(raw::AbstractString)
    section=match(r"\"cross_bus_candidates\"\s*:\s*\[(.*?)\]\s*,\s*\"domain\""s,raw)
    section===nothing && error("candidate cross_bus_candidates array is missing")
    rows=NamedTuple[]
    for item in eachmatch(r"\{([^{}]*)\}",section.captures[1])
        o=item.captures[1]
        push!(rows,(bus=round(Int,json_number(o,"bus")),rho=json_number(o,"rho_star"),
            beta=json_number(o,"pll_scale"),kp=json_number(o,"Kp_star"),ki=json_number(o,"Ki_star"),
            alpha=json_number(o,"alpha"),pole=complex(json_number(o,"rightmost_pole_real"),
                json_number(o,"rightmost_pole_imag")),all_modes_pass=json_bool(o,"all_modes_pass")))
    end
    isempty(rows) && error("frozen candidate contains no bus candidates")
    return rows
end

"""Validate frozen candidates and evaluate a post-freeze rho/PLL grid.

The candidate and sidecar are checked before PowerDynamics is imported. The
validation grid follows the preregistered physical PLL coordinate
`Kp=beta*Kp0`, `Ki=beta^2*Ki0`; it evaluates every point on a fixed 11-by-3
rho/beta grid for each bus, using the PowerDynamics current-sharing model.
"""
function run_validation(root)
    gate=candidate_gate(root)
    gate.allowed || return (status="REFUSED_NO_FROZEN_CANDIDATE",reason=gate.reason,
        powerdynamics_imported=false,powerdynamics_calls=0)

    # The import and every network call occur only after hash/status verification.
    Core.eval(@__MODULE__, :(using PowerDynamics))
    Core.eval(@__MODULE__, :(using NetworkDynamics))
    Base.include(@__MODULE__,joinpath(root,"src","pd39","PD39.jl"))
    pd39=getfield(@__MODULE__,:PD39)
    pdmodel=getfield(pd39,:PD39Model)

    rows=NamedTuple[]
    calls=0
    sigma=json_number(gate.raw,"sigma_req")
    kp0=5.0*2pi
    ki0=kp0^2/4
    frozen=candidate_rows(gate.raw)
    endpoint_records=Dict{Tuple{Int,Float64},NamedTuple}()
    endpoint_spectra=Dict{Tuple{Int,Float64},Vector{ComplexF64}}()
    for c in frozen
        c.rho==1.0 || error("endpoint validator accepts only the frozen saturated-rho candidates")
        c.all_modes_pass || error("analytic all-mode gate is false for bus $(c.bus)")
        isapprox(c.kp,kp0*c.beta;rtol=1e-12,atol=1e-12) || error("candidate Kp is inconsistent with PLL scale")
        isapprox(c.ki,ki0*c.beta^2;rtol=1e-12,atol=1e-12) || error("candidate Ki is inconsistent with PLL scale")

        nw0=Base.invokelatest(pdmodel.baseline_network)
        template=Base.invokelatest(pdmodel.simple_gfldc_template;pll_scale=c.beta)
        nw=Base.invokelatest(pdmodel.replace_bus,nw0,c.bus;template)
        eq=Base.invokelatest(pd39.initialize_equilibrium,nw;sparse=false,check=:error)
        calls+=1
        λ=ComplexF64.(Base.invokelatest(NetworkDynamics.jacobian_eigenvals,eq.state))
        endpoint_spectra[(c.bus,c.beta)]=λ
        finite=all(isfinite,real.(λ))&&all(isfinite,imag.(λ))
        nongauge=filter(z->abs(z)>1e-8,λ)
        isempty(nongauge)&&error("PowerDynamics returned no non-gauge poles at bus $(c.bus)")
        pdpole=nongauge[argmax(real.(nongauge))]
        pdalpha=real(pdpole)
        nearest=nongauge[argmin(abs.(nongauge .- c.pole))]
        pole_error=abs(nearest-c.pole)
        alpha_error=abs(pdalpha-c.alpha)
        audit=Base.invokelatest(pd39.stability_audit,eq.state)
        pass=eq.powerflow_finite&&eq.state_finite&&eq.fixed_point&&finite&&
            pdalpha<=-sigma+1e-10&&pole_error<=1e-5&&alpha_error<=1e-6
        record=(bus=c.bus,rho=c.rho,beta=c.beta,Kp=c.kp,Ki=c.ki,
            analytic_rightmost_real=real(c.pole),analytic_rightmost_imag=imag(c.pole),
            PD_rightmost_real=real(pdpole),PD_rightmost_imag=imag(pdpole),
            analytic_alpha=c.alpha,PD_alpha=pdalpha,alpha_abs_error=alpha_error,
            nearest_pole_abs_error=pole_error,PD_margin=audit.dynamic_margin,
            equilibrium_pass=eq.powerflow_finite&&eq.state_finite&&eq.fixed_point,
            non_gauge_poles=length(nongauge),all_modes_pass=all(real.(nongauge).<=-sigma+1e-10),
            status=pass ? "PASS" : "FAIL")
        push!(rows,record)
        endpoint_records[(c.bus,c.beta)]=record
    end
    all(r.status=="PASS" for r in rows) || error("independent PowerDynamics candidate validation failed")

    # One baseline PF/equilibrium is shared because the current-share model
    # preserves the frozen bus injections. Every interior gain/share point is
    # initialized independently from a copy of these PF data.
    rho_grid=collect(0.0:0.1:1.0)
    beta_grid=[0.9,1.0,1.1]
    baseline_nw=Base.invokelatest(pdmodel.baseline_network)
    baseline_eq=Base.invokelatest(pd39.initialize_equilibrium,baseline_nw;sparse=false,check=:error)
    calls+=1
    baseline_lambda=ComplexF64.(Base.invokelatest(NetworkDynamics.jacobian_eigenvals,baseline_eq.state))
    baseline_finite=all(isfinite,real.(baseline_lambda))&&all(isfinite,imag.(baseline_lambda))
    baseline_nongauge=filter(z->abs(z)>1e-8,baseline_lambda)
    isempty(baseline_nongauge)&&error("PowerDynamics baseline has no non-gauge poles")
    baseline_pole=baseline_nongauge[argmax(real.(baseline_nongauge))]
    baseline_fixed=baseline_eq.powerflow_finite&&baseline_eq.state_finite&&baseline_eq.fixed_point&&baseline_finite
    baseline_alpha=real(baseline_pole)
    baseline_pass=baseline_fixed&&baseline_alpha<=-sigma+1e-10

    grid_rows=NamedTuple[]
    full_spectra_rows=NamedTuple[]
    for bus in sort(unique(c.bus for c in frozen))
        # Compile the 121-state interior topology once; change only the frozen
        # share and PLL defaults at subsequent grid points.
        interior_template=Base.invokelatest(pdmodel.weighted_replacement_network,baseline_nw,bus,0.5,1.0)
        for rho in rho_grid, beta in beta_grid
            kp=kp0*beta
            ki=ki0*beta^2
            eqpass=false; finite=false; alpha=NaN; pole=complex(NaN,NaN); nmodes=0
            mode_values=ComplexF64[]
            status="EQUILIBRIUM_OR_SPECTRUM_FAILURE"; margin=NaN; reused=false
            try
                if rho==0.0
                    eqpass=baseline_fixed; finite=baseline_finite; alpha=baseline_alpha
                    mode_values=baseline_lambda
                    pole=baseline_pole; nmodes=length(baseline_nongauge)
                    margin=-alpha; status=baseline_pass ? "PASS" : "SPECTRAL_VIOLATION"
                    reused=true
                elseif rho==1.0 && beta==1.0 && haskey(endpoint_records,(bus,beta))
                    e=endpoint_records[(bus,beta)]
                    eqpass=e.equilibrium_pass; finite=true; alpha=e.PD_alpha
                    mode_values=endpoint_spectra[(bus,beta)]
                    pole=complex(e.PD_rightmost_real,e.PD_rightmost_imag)
                    nmodes=e.non_gauge_poles; margin=e.PD_margin
                    status=e.all_modes_pass ? "PASS" : "SPECTRAL_VIOLATION"
                    reused=true
                else
                    nw = if rho==1.0
                        template=Base.invokelatest(pdmodel.simple_gfldc_template;pll_scale=beta)
                        Base.invokelatest(pdmodel.replace_bus,baseline_nw,bus;template)
                    else
                        Base.invokelatest(pdmodel.set_weighted_replacement_parameters,interior_template,bus,rho,beta)
                    end
                    eq=Base.invokelatest(pd39.initialize_equilibrium,nw;pfs=deepcopy(baseline_eq.pfs),sparse=false,check=:error)
                    calls+=1
                    λ=ComplexF64.(Base.invokelatest(NetworkDynamics.jacobian_eigenvals,eq.state))
                    mode_values=λ
                    finite=all(isfinite,real.(λ))&&all(isfinite,imag.(λ))
                    nongauge=filter(z->abs(z)>1e-8,λ)
                    nmodes=length(nongauge)
                    isempty(nongauge)&&error("no non-gauge poles")
                    pole=nongauge[argmax(real.(nongauge))]
                    alpha=real(pole); margin=-alpha
                    eqpass=eq.powerflow_finite&&eq.state_finite&&eq.fixed_point
                    status=eqpass&&finite ? (alpha<=-sigma+1e-10 ? "PASS" : "SPECTRAL_VIOLATION") : "EQUILIBRIUM_OR_SPECTRUM_FAILURE"
                end
            catch err
                status="EQUILIBRIUM_OR_SPECTRUM_FAILURE"
                alpha=NaN; pole=complex(NaN,NaN); margin=NaN
                eqpass=false; finite=false; nmodes=0
                mode_values=ComplexF64[]
                @warn "PowerDynamics ExpD grid cell failed" bus rho beta exception=(err,catch_backtrace())
            end
            push!(grid_rows,(bus=bus,rho=rho,beta=beta,Kp=kp,Ki=ki,
                rightmost_real=real(pole),rightmost_imag=imag(pole),alpha=alpha,
                margin=margin,equilibrium_pass=eqpass,spectrum_finite=finite,
                non_gauge_poles=nmodes,within_sigma=status=="PASS",reused_evaluation=reused,
                status=status))
            for (mode,λj) in enumerate(mode_values)
                gauge=abs(λj)<=1e-8
                push!(full_spectra_rows,(bus=bus,rho=rho,beta=beta,Kp=kp,Ki=ki,mode=mode,
                    lambda_real=real(λj),lambda_imag=imag(λj),is_gauge=gauge,
                    mode_pass=gauge||real(λj)<=-sigma+1e-10,
                    source=(rho==0.0 ? "PD_baseline_reused" :
                        (rho==1.0&&beta==1.0 ? "PD_frozen_candidate_reused" : "PD_grid_linearization"))))
            end
        end
    end
    grid_complete=length(grid_rows)==length(unique(c.bus for c in frozen))*length(rho_grid)*length(beta_grid) &&
        all(r.status!="EQUILIBRIUM_OR_SPECTRUM_FAILURE" for r in grid_rows)

    table=DataFrame(rows)
    grid_table=DataFrame(grid_rows)
    full_spectra_table=DataFrame(full_spectra_rows)
    outdir=joinpath(root,"reports","experiment_D","tables")
    mkpath(outdir)
    CSV.write(joinpath(outdir,"TABLE_D11_PD_validation.csv"),table)
    CSV.write(joinpath(outdir,"TABLE_D17_PD_rho_PLL_grid.csv"),grid_table)
    CSV.write(joinpath(outdir,"TABLE_D19_PD_grid_full_spectra.csv"),full_spectra_table)
    open(joinpath(root,"reports","experiment_D","PD_VALIDATION_EXP_D.json"),"w") do io
        print(io,"{\"status\":\"PASS\",\"powerdynamics_calls\":",calls,
            ",\"candidate_sha256\":\"",gate.sha256,"\",\"grid_status\":\"",
            grid_complete ? "COMPLETE" : "INCOMPLETE","\",\"grid_rows\":",length(grid_rows),
            ",\"grid_feasible_rows\":",count(r->r.within_sigma,grid_rows),
            ",\"grid_full_spectrum_rows\":",length(full_spectra_rows),"}\n")
    end
    return (status="PASS",reason="all frozen endpoints match a fresh PowerDynamics linearization; post-freeze rho/PLL grid evaluated",
        candidate_sha256=gate.sha256,powerdynamics_imported=true,powerdynamics_calls=calls,
        table=table,grid_table=grid_table,full_spectra_table=full_spectra_table,
        grid_status=grid_complete ? "COMPLETE" : "INCOMPLETE")
end

end
