using LinearAlgebra,CSV,DataFrames,TOML,Statistics
include(joinpath(@__DIR__,"DDE_NEV.jl"))
const N=DDENEV;const D=N.D;const R=D.R;const BASE=joinpath(@__DIR__,"baseline_reproduction")
BLAS.set_num_threads(1);ctx=R.N.design_context(R.ROOT)
rho=fill(.875,10);kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
L=D.linearization(ctx,rho,kp,ki)
frozen=CSV.read(joinpath(@__DIR__,"FROZEN_DELAY_PATTERNS.csv"),DataFrame)
cases=[("uniform_10ms",fill(10.0,10)),
       ("min_commutator",parse.(Float64,split(only(frozen.tau_ms[frozen.delay_pattern_id.=="min_commutator"]),";"))),
       ("max_commutator",parse.(Float64,split(only(frozen.tau_ms[frozen.delay_pattern_id.=="max_commutator"]),";")))]
pole_rows=NamedTuple[];count_rows=NamedTuple[]
for (case,tms) in cases
    tau=tms./1000;println("DDE_ROOTS_START ",case," delay_ms=",tms);flush(stdout)
    tracked=N.track_roots(L,tau;real_cut=-.25,continuation_steps=4,max_roots=2,maxiter=12)
    inside=[r for r in tracked if r.converged && real(r.s)>-.25]
    match_count=false
    push!(count_rows,(;delay_pattern_id=case,mean_tau_ms=mean(tms),gamma=-.25,
      norm_derived_radius=missing,contour_root_count=missing,tracked_roots_inside=length(inside),
      tracked_root_count=length(tracked),converged_root_count=count(r->r.converged,tracked),
      max_phase_increment=missing,contour_points=missing,
      minimum_logabs_small_determinant=missing,contour_status="NOT_COMPLETED_CONSERVATIVE_BOUND_COST",
      root_coverage_match=match_count,status="BLOCKED_EXACT_DDE_SPECTRUM"))
    for r in tracked
        z=r.s
        push!(pole_rows,(;case_id="baseline",delay_pattern_id=case,tau_vector=join(tms,";"),pole_id=r.root_id,
          real=real(z),imag=imag(z),frequency_hz=abs(imag(z))/(2pi),damping_ratio=-real(z)/abs(z),
          residual=r.residual,solver="exact_characteristic_bordered_Newton",root_tracking_id=r.root_id,
          converged=r.converged,continuation_steps=4,contour_count=missing,root_coverage_match=match_count))
    end
    println("DDE_ROOTS_DONE ",case," tracked=",length(tracked)," converged=",count(r->r.converged,tracked),
      " inside=",length(inside)," status=BLOCKED_EXACT_DDE_SPECTRUM");flush(stdout)
end
CSV.write(joinpath(BASE,"TABLE_D01_CONTOUR_COUNTS.csv"),DataFrame(count_rows))
CSV.write(joinpath(BASE,"TABLE_D01_TRACKED_DDE_POLES.csv"),DataFrame(pole_rows))
zero=CSV.read(joinpath(BASE,"TABLE_D01_EXACT_DDE_POLES.csv"),DataFrame)
tracked=DataFrame(pole_rows)
zero.pole_id=string.(zero.pole_id)
append!(zero,tracked;cols=:union)
CSV.write(joinpath(BASE,"TABLE_D01_EXACT_DDE_POLES.csv"),zero)
println(DataFrame(count_rows))
