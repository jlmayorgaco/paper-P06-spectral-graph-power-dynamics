using CSV, DataFrames, DelimitedFiles, LinearAlgebra, NetworkDynamics, PowerDynamics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_N","identity_probe_bus38_half")
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))

function save_matrix(prefix,sys)
    n=size(sys.A,1)
    m=sys.M isa UniformScaling ? Matrix{Float64}(I,n,n) : Matrix(sys.M)
    for (name,val) in (("M",m),("A",sys.A),("B",sys.B),("C",sys.C),("D",sys.D))
        val===nothing && continue
        writedlm(prefix*"_"*name*".csv",Matrix(val),',')
    end
end

function main()
    mkpath(OUT)
    println("N_IDENTITY_BASELINE");flush(stdout)
    base=PDReferenceN.frozen_baseline()
    ctx=PDExactDesignN.design_context(ROOT)
    rho=zeros(10);rho[38-29]=0.5
    kp=fill(PDExactDesignN.K0P,10)
    ki=fill(PDExactDesignN.K0I,10)
    println("N_IDENTITY_PD_BUILD");flush(stdout)
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    s=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    res=PDReferenceN.residual_audit(nw,s)
    println("N_IDENTITY_LINEARIZE residual=",res.maximum);flush(stdout)
    sys=linearize_network(s)
    save_matrix(joinpath(OUT,"PD"),sys)
    mass=sys.M isa UniformScaling ? ones(size(sys.A,1)) : diag(sys.M)
    CSV.write(joinpath(OUT,"PD_state_map.csv"),DataFrame(
        index=1:length(sys.sym),state_name=string.(sys.sym),
        differential=mass.==1.0))
    λpd=jacobian_eigenvals(sys)
    CSV.write(joinpath(OUT,"PD_poles.csv"),DataFrame(real=real.(λpd),imag=imag.(λpd)))
    buslin=linearize_component(s,VIndex(38))
    save_matrix(joinpath(OUT,"PD_bus38"),buslin)
    println("N_IDENTITY_ANALYTIC");flush(stdout)
    ans=PDExactDesignN.spectrum(ctx,rho,kp,ki)
    m=ans.model
    for (name,val) in (("Ared",m.Ared),("A_local",m.A),("B_port",m.B),
                       ("C_port",m.C),("D_port",m.D))
        writedlm(joinpath(OUT,"AN_"*name*".csv"),Matrix(val),',')
    end
    CSV.write(joinpath(OUT,"AN_state_map.csv"),m.state_inventory)
    CSV.write(joinpath(OUT,"AN_poles.csv"),DataFrame(real=real.(ans.lambda),imag=imag.(ans.lambda)))
    println("N_IDENTITY_DONE PD_count=",length(λpd)," AN_count=",length(ans.lambda),
        " PD_alpha_raw=",maximum(real.(λpd))," AN_alpha=",ans.alpha)
end

main()
