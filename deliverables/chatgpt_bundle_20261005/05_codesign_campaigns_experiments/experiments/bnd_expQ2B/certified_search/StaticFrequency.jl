module StaticFrequency
using LinearAlgebra
import ..LocalOracle as L
function build()
    ctx=L.CTX;gen=collect(59:78);rest=collect(1:58);Y=ctx.net.y_static
    kron=Y[gen,gen]-Y[gen,rest]*(Y[rest,rest]\Y[rest,gen])
    dummy=L.N.descriptor(ctx,fill(.5,10));bd=L.LinearSecurity.load_input_vector(ctx,dummy,16)
    b=bd[gen]-Y[gen,rest]*(Y[rest,rest]\bd[rest])
    A0=zeros(21,21);A0[1:20,1:20]=kron;delta=Matrix{Float64}[]
    for i in 1:10
        bus=i+29;rows=2i-1:2i;op=L.N.trim_gfl(ctx,bus)
        gf=L.N.gfl_jacobians(op,L.N.K0P,L.N.K0I);sg=ctx.net.sg[bus].J
        gg=zeros(9);gg[3]=1;gg[6]=op.x[7];gg[7]=-op.x[6]
        gs=zeros(size(sg.A,1));gs[end]=1
        yg=gf.D-gf.C*(gf.A\gf.B);fg=gf.C*(gf.A\gg)*2pi
        ys=sg.D-sg.C*(sg.A\sg.B);fs=sg.C*(sg.A\gs)*2pi
        A0[rows,rows]+=yg;A0[rows,21]=fg
        d=zeros(21,21);d[rows,rows]=ys-yg;d[rows,21]=fs-fg;push!(delta,d)
        u=ctx.net.voltage[bus];A0[21,rows]=[-imag(u),real(u)]
    end
    rhs=vcat(-100b,0.)
    # Fixed row equilibration preserves both determinant ratios and parameter rank.
    scales=[max(maximum(abs.(A0[r,:])),maximum(maximum(abs.(d[r,:])) for d in delta),1e-12) for r in 1:21]
    A0=A0./scales;delta=[d./scales for d in delta];rhs=rhs./scales
    (;A0,delta,rhs)
end
matrix(m,e)=m.A0+sum(e[i]*m.delta[i] for i in 1:10)
frequency(m,e)=(matrix(m,e)\m.rhs)[end]
function determinants(m,e)
    A=matrix(m,e);D=det(A);A[:,end]=m.rhs
    (;D,N=det(A))
end
end
