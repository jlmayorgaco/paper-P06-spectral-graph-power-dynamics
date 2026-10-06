include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, Random, Test, TOML
const R=ReducedDAE
BLAS.set_num_threads(1);ctx=R.N.design_context(R.ROOT);rng=MersenneTwister(61497)
maxerr=0.;legacy_max=0.
@testset "Converter plus AC filter energy balance" begin
    for bus in 30:39,trial in 1:20
        op=R.N.trim_gfl(ctx,bus);p=op.pars;x=copy(op.x)
        x.+=1e-5randn(rng,9);u=op.u+.001randn(rng,2)
        phys=R.gfl_rhs(x,u,p,R.N.K0P,R.N.K0I,:physical_supply)
        legacy=R.gfl_rhs(x,u,p,R.N.K0P,R.N.K0I,:legacy)
        # A single common system base; port_scale would multiply both stored energies
        # and terminal/DC powers for the aggregate, preserving the identity.
        energy_dot=p.Cdc*x[9]*phys[9]+p.Xf/p.omega_base*(x[6]*phys[6]+x[7]*phys[7])
        supply=p.Pdc-(u[1]*x[7]+u[2]*x[6])-p.Rf*(x[6]^2+x[7]^2)
        err=abs(energy_dot-supply);global maxerr=max(maxerr,err)
        @test err<1e-10
        wrong=p.Cdc*x[9]*legacy[9]+p.Xf/p.omega_base*(x[6]*legacy[6]+x[7]*legacy[7])
        global legacy_max=max(legacy_max,abs(wrong-supply))
    end
end
rho=fill(.8,10);kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
ml=R.model(ctx,rho,kp,ki);mp=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
Jl=R.derivatives(ml.x0,ml).Fx;Jp=R.derivatives(mp.x0,mp).Fx
S=ones(length(mp.x0));for ix in mp.gfidx;S[ix[9]]=-1.;end
similarity=norm(Jp-Diagonal(S)*Jl*Diagonal(S),Inf)/norm(Jl,Inf)
@test similarity<1e-12
trim=zeros(length(mp.x0));R.rhs!(trim,mp.x0,mp,0.);@test norm(trim,Inf)<1e-8
out=Dict("physical_energy_balance_max_error_pu"=>maxerr,"legacy_physical_energy_balance_max_error_pu"=>legacy_max,
    "sample_count"=>200,"nominal_Jacobian_similarity_error"=>similarity,"physical_trim_residual"=>norm(trim,Inf),
    "interpretation"=>"With positive converter-to-AC power and positive DC supply, physical balance is d(E_dc+E_filter)/dt=Pdc-Pbus-Rf*I^2. Both the inherited DC balance and voltage-control errors had opposite signs. Corrected physical convention is a separate model variant. At trim the two linearizations are similar under delta-vdc sign reversal; nonlinear equivalence is not asserted.")
mkpath(R.OUT);open(joinpath(R.OUT,"dc_physics_audit.toml"),"w") do io;TOML.print(io,out);end
println("DC_PHYSICS_AUDIT ",out);flush(stdout)
