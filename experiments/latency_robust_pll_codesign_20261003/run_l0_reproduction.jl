using CSV, DataFrames, LinearAlgebra, TOML

module Roots
include(joinpath(@__DIR__, "..", "optimal_latency_margin_20261002", "run_fast_root_locus.jl"))
end

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PREV = joinpath(ROOT, "experiments", "optimal_latency_margin_20261002")
const OUT = @__DIR__
BLAS.set_num_threads(1)

function classify(L, tau_ms)
    q = Roots.MegaOracle.trace_count(L, fill(tau_ms / 1000, length(L.Ai)); gamma=-0.05, step=10.0)
    good = abs(real(q.estimate)-q.nearest)<0.01 && abs(imag(q.estimate))<0.01 &&
        abs(q.phase_estimate-q.nearest)<0.01 && q.largest_phase_step<pi/2 &&
        q.quadrature_error_estimate<0.01
    (; status=good ? (q.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
       count=good ? q.nearest : missing,
       trace_count_real=real(q.estimate), phase_count=q.phase_estimate,
       phase_step=q.largest_phase_step, quadrature_error=q.quadrature_error_estimate,
       min_sigma=q.min_sigma_small, evaluations=q.evaluations)
end

function main()
    ctx=Roots.MegaOracle.DC.R.N.design_context(ROOT)
    prevZ=CSV.read(joinpath(PREV,"T02_PRECISE_CROSSINGS.csv"),DataFrame)
    prevN=CSV.read(joinpath(PREV,"T04_PRECISE_CROSSINGS.csv"),DataFrame)
    root_seeds=Dict("Z"=>prevZ[prevZ.design_id.=="best_zero_delay_88455",:][1,:],
                    "N"=>prevN[prevN.design_id.=="rho_04",:][1,:])
    rows=NamedTuple[]; countrows=NamedTuple[]; modevectors=Dict{String,Vector{ComplexF64}}()
    designs=(("Z",joinpath(OUT,"designs","Z_zero_delay_tuned.toml")),
             ("N",joinpath(OUT,"designs","N_nominal.toml")))
    for (id,path) in designs
        d=TOML.parsefile(path);rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
        L=Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
        model=Roots.reduced_model(L)
        spectrum=eigvals(L.A)
        iz=argmax(real.(spectrum));alpha0=real(spectrum[iz]);f0=abs(imag(spectrum[iz]))/(2pi)
        seed=root_seeds[id]
        s0=id=="Z" ? ComplexF64(seed.critical_root_real+im*seed.critical_root_imag) :
                       ComplexF64(seed.critical_root_real+im*(2pi*seed.critical_frequency_hz))
        function at(tau_ms)
            r=Roots.refine(model,s0,fill(tau_ms/1000,length(L.Ai));tol=1e-11)
            r.converged || error("root refinement failed $id $tau_ms ms")
            r
        end
        lo,hi= id=="Z" ? (37.36,37.42) : (39.35,39.42)
        rlo,rhi=at(lo),at(hi)
        println("L0_INITIAL_ROOT_BRACKET ",id," lo=",lo," real=",real(rlo.s),
            " hi=",hi," real=",real(rhi.s));flush(stdout)
        real(rlo.s)<-0.05<real(rhi.s) || error("root not bracketed for $id")
        while hi-lo>0.002
            mid=(lo+hi)/2
            if real(at(mid).s)<-0.05;lo=mid;else;hi=mid;end
        end
        tau=(lo+hi)/2;r=at(tau)
        D=Roots.MegaOracle.DC.delta_matrix(L,r.s,fill(tau/1000,length(L.Ai)))
        S=svd(D);v=ComplexF64.(L.Q*S.V[:,end]);v/=norm(v);modevectors[id]=v
        samples=unique([20.0,30.0,35.0,tau-0.02,tau+0.02,40.0])
        below_status="";above_status="";above_count=missing
        for x in samples
            println("L0_COUNT_START ",id," tau_ms=",x);flush(stdout)
            try
                c=classify(L,x)
                push!(countrows,merge((;design_id=id,tau_ms=x),c))
                if abs(x-(tau-0.02))<1e-8;below_status=c.status;end
                if abs(x-(tau+0.02))<1e-8;above_status=c.status;above_count=c.count;end
            catch err
                push!(countrows,(;design_id=id,tau_ms=x,status="INDETERMINATE",count=missing,
                    trace_count_real=NaN,phase_count=NaN,phase_step=NaN,quadrature_error=NaN,
                    min_sigma=NaN,evaluations=missing))
                println("L0_COUNT_ERROR ",id," ",x," ",sprint(showerror,err));flush(stdout)
            end
            CSV.write(joinpath(OUT,"L0_ROOT_COUNTS.csv"),DataFrame(countrows))
        end
        row=(;design_id=id,rho_i=join(rho,";"),Kp_i=join(kp,";"),Ki_i=join(ki,";"),
            alpha_zero_delay_s_inv=alpha0,critical_zero_delay_frequency_hz=f0,
            tau_crit_local_ms=tau,tau_crit_local_bracket_ms=hi-lo,
            critical_delay_frequency_hz=imag(r.s)/(2pi),critical_root_real=real(r.s),
            critical_root_imag=imag(r.s),root_residual=r.sv,
            below_status,above_status,above_margin_violating_roots=above_count,
            critical_family_id=id=="Z" ? "Z_fast5" : "N_fast3",
            zero_delay_event_status="PENDING_LOCAL_RECOMPUTATION",
            status="LOCAL_ROOT_RECOMPUTED_EXACT_CHARACTERISTIC;FULL_CONTOUR_BOUNDARY_CHECK")
        push!(rows,row);CSV.write(joinpath(OUT,"L0_REPRODUCTION.csv"),DataFrame(rows))
        println("L0_RESULT ",row);flush(stdout)
    end
    overlap=abs(dot(modevectors["Z"],modevectors["N"]))^2
    CSV.write(joinpath(OUT,"L0_MODAL_OVERLAP.csv"),DataFrame([(;design_left="Z",design_right="N",
        full_physical_right_vector_MAC=overlap,status="RIGHT_VECTOR_MAC_NOT_BIORTHOGONAL_PARTICIPATION")]))
    println("L0_MODAL_MAC ",overlap)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
