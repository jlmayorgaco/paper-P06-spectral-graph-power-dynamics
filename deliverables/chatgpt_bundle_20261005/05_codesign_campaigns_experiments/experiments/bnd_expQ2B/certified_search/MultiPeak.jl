# Included inside LocalOracle: semi-infinite constraint exchange, keeping
# distinct temporal maxima instead of differentiating their nonsmooth max.
function peak_surfaces(v,kind;threshold=.35)
    grid=sort(unique(vcat(0.,10. .^range(-6,log10(.02),length=32),collect(.02:.04:60.),.5-1e-10,.5,.5+1e-10,1-1e-10,1.,1+1e-10)))
    HT=reduce(hcat,[coeff_t(v.lam,t,kind) for t in grid])
    ddt=kind==:F ? Float64.(grid.<.5)./.5 : (Float64.(grid.<.5).-Float64.((grid.>=.5).&(grid.<1)))./.25
    dvalues=100real.(v.R*HT+v.D*ddt')
    late=findall(t->t>=(kind==:F ? .5 : 1.),grid)
    tail=reduce(hcat,[tailcoeff(v.lam,grid[k],kind) for k in late])
    dvalues[:,late]=100real.(v.R*(tail.*v.lam))
    result=NamedTuple[]
    for ch in 1:10
        times=[0.,.5-1e-10,.5,.5+1e-10,1-1e-10,1.,1+1e-10]
        dt=dvalues[ch,:]
        for j in 1:length(grid)-1
            (grid[j]<.5<=grid[j+1] || grid[j]<1<=grid[j+1]) && continue
            signbit(dt[j])==signbit(dt[j+1]) && continue
            root=bisect(t->output_dt(v,t,kind,ch),grid[j],grid[j+1])
            output(v,root,kind,ch)*dt[j]>0 && push!(times,root)
        end
        for t in times
            f=output(v,t,kind,ch)
            abs(f)<threshold && continue
            push!(result,(;value=abs(f),time=t,sign=sign(f),channel=ch,kind))
        end
    end
    sort!(result;by=p->(p.channel,p.time))
end
function surface_gradient(v,p)
    lam=v.lam;V=v.V;Vi=v.Vi;Z=v.Z;CV=v.CV
    late=(p.kind==:F && p.time>=.5)||(p.kind==:R && p.time>=1.)
    fn=late ? tailcoeff : coeff
    fs=fn(lam,p.time,p.kind);dfs=fn(lam,p.time,p.kind;derivative=true)
    dd=[abs(lam[i]-lam[j])<1e-7*max(1.,abs(lam[i])) ? (dfs[i]+dfs[j])/2 :
        (fs[i]-fs[j])/(lam[i]-lam[j]) for i in eachindex(lam),j in eachindex(lam)]
    ell=(reshape(CV[p.channel,:],:,1).*dd).*reshape(Z,1,:)
    weight=transpose(Vi)*ell*transpose(V)
    d,jump=late ? (0.,0.) : coefficients(p.time,p.kind)
    ab=v.A\v.B
    [begin
        val=sum((z.C[p.channel,:]'*V)[:].*fs.*Z)+sum(CV[p.channel,:].*fs.*(Vi*z.B))+
            sum(weight.*z.A)+z.D[p.channel]*d+z.jump[p.channel]*jump
        dc_z=late && p.kind==:F ? 100real((-z.C*ab+v.C*(v.A\(z.A*ab-z.B))+z.D)[p.channel]) : 0.
        p.sign*(100real(val)+dc_z)/.5
    end for z in v.ds]
end
function multi_evaluate(a,x;derivatives=false)
    v=evaluate(a,x;derivatives)
    ps=vcat(peak_surfaces(v,:F),peak_surfaces(v,:R))
    # Same steady frequency at connected buses; keep one exact independent row.
    rows=[1,2,3]
    gs=vcat(v.g[rows],[p.value/.5-1 for p in ps])
    names=vcat(["modal","beta","steady"],["$(p.kind)_bus$(p.channel+29)_t$(round(p.time;digits=5))" for p in ps])
    derivatives || return merge(v,(;g=gs,names,surfaces=ps))
    J=vcat(v.Jac[rows,:],reduce(vcat,[surface_gradient(v,p)' for p in ps];init=zeros(0,length(x))))
    merge(v,(;g=gs,Jac=J,names,surfaces=ps))
end
function matched_jacobian(vnew,vold)
    J=copy(vold.Jac);J[1:3,:]=vnew.Jac[1:3,:]
    for (i,p) in enumerate(vold.surfaces)
        choices=findall(q->q.channel==p.channel && q.kind==p.kind,vnew.surfaces)
        if isempty(choices)
            J[i+3,:]=surface_gradient(vnew,p)
        else
            k=choices[argmin([abs(vnew.surfaces[k].time-p.time) for k in choices])]
            J[i+3,:]=vnew.Jac[k+3,:]
        end
    end
    J
end
