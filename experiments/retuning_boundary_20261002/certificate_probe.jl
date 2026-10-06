include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
using .ReducedDAE,LinearAlgebra,CSV,DataFrames,TOML
const R=ReducedDAE
const OUT=joinpath(R.ROOT,"reports","retuning_boundary_20261002")
BLAS.set_num_threads(1)
say(x...)=(println(x...);flush(stdout))
const CTX=R.N.design_context(R.ROOT)

function pencil(rho,kp,ki)
    m=R.model(CTX,rho,kp,ki;delta=0.,dc_convention=:physical_supply)
    de=R.derivatives(m.x0,m);q,_=R.rotation_generator(m.x0,m)
    Q=nullspace(reshape(q/norm(q),1,:));nx=length(m.x0)
    U=zeros(nx,20);C=zeros(20,nx)
    for i in 1:10
        ix=m.gfidx[i];theta=m.x0[ix[3]];s,c=sincos(theta)
        measurement=-s*de.vx[2i-1,:]+c*de.vx[2i,:]
        measurement[ix[3]]-=c*de.v[2i-1]+s*de.v[2i]
        C[i,:]=measurement;C[10+i,:]=measurement
        U[ix[4],i]=1/m.gp[i].pll_tau
        U[ix[5],10+i]=1.
    end
    (;m,A=Q'*de.Fx*Q,U=Q'*U,C=C*Q,Q)
end

function scan()
    rows=NamedTuple[]
    for rho in (.875,.90,.925,.95,.975,.99)
        p=pencil(fill(rho,10),fill(R.N.K0P,10),fill(R.N.K0I,10))
        vals=eigvals(p.A);j=argmax(real.(vals))
        push!(rows,(;rho,alpha=real(vals[j]),imag=imag(vals[j])))
        say("SCAN ",last(rows))
    end
    CSV.write(joinpath(OUT,"preflight_spectrum.csv"),DataFrame(rows))
    for rho in (.925,.95,.975,.99),width in (.1,.5,.9)
        center=width==.9 ? 2.125 : 1.
        pp=pencil(fill(rho,10),fill(center*R.N.K0P,10),fill(center*R.N.K0I,10))
        ee=eigen(pp.A);j=argmax(real.(ee.values));lambda=ee.values[j]
        real(lambda)>-.05 || continue
        radius=min((real(lambda)+.05)*.7,.25)
        radii=vcat(fill(width*center*R.N.K0P,10),fill(width*center*R.N.K0I,10))
        B=ee.vectors\(pp.U*Diagonal(radii));C=pp.C*ee.vectors
        maxnorm=0.;maxradius=0.
        for theta in range(0,2pi;length=129)[1:end-1]
            s=lambda+radius*cis(theta)
            M=C*(B./(s.-ee.values))
            maxnorm=max(maxnorm,opnorm(M))
            maxradius=max(maxradius,maximum(abs.(eigvals(M))))
        end
        say("CERT_PROBE rho=",rho," width=",width," center=",center,
            " lambda=",lambda," radius=",radius," norm=",maxnorm," spectral=",maxradius)
    end
end
scan()
