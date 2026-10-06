include("setup.jl")
using Sockets
const S=collect(30:39)
const M0=N.descriptor(CTX,ORIGINAL["rho"],ORIGINAL["Kp"],ORIGINAL["Ki"])
const G=N.gauge_vector(M0)
const Q=nullspace(reshape(G/norm(G),1,:))
const ANG=let W=zeros(10,78)
    for bus in 30:39
        v=CTX.net.voltage[bus];W[bus-29,2bus-1]=-imag(v)/abs2(v)
        W[bus-29,2bus]=real(v)/abs2(v)
    end;W
end
const LOADS=hcat([LinearSecurity.load_input_vector(CTX,M0,b) for b in (16,8,29)]...)

function dynamic_jacobian(y)
  v=oracle(y);A=v.A
  er=eigen(A);j=argmax(real.(er.values));lambda=er.values[j];r=er.vectors[:,j]
  el=eigen(A');l=el.vectors[:,argmin(abs.(el.values.-conj(lambda)))];den=dot(l,r)
  F=svd(im*v.omega*I-A-.05I);u=F.U[:,end];w=F.V[:,end]
  ep,kp,ki=CoreDesign.decode(CTX,S,y);m=N.descriptor(CTX,1 .-ep,kp,ki)
  gv=N.gauge_vector(m);qb=nullspace(reshape(gv/norm(gv),1,:))
  function getA(x)
    e,p,i=CoreDesign.decode(CTX,S,x);mm=N.descriptor(CTX,1 .-e,p,i)
    qb'*mm.Ared*qb
  end
  J=zeros(2,30)
  for col in 1:30
    col<=10 && y[col]==0 && continue
    h=col<=10 ? min(2e-6,y[col]/4,(1-y[col])/4) : 2e-6
    xp=copy(y);xm=copy(y);xp[col]=min(1-1e-12,y[col]+h);xm[col]=max(col<=10 ? 1e-10 : 0.0,y[col]-h)
    Az=(getA(xp)-getA(xm))/(xp[col]-xm[col])
    J[1,col]=real(dot(l,Az*r)/den)/.05
    if 1-v.beta/(1.02CoreDesign.BETA_REQ)>-5
      J[2,col]=-real(dot(u,-Az*w))/(1.02CoreDesign.BETA_REQ)
    end
  end
  J
end

function oracle(y;dt=0.05,horizon=40.0,extra_frequencies=Float64[])
    ep,kp,ki=CoreDesign.decode(CTX,S,y);rho=1 .-ep
    m=N.descriptor(CTX,rho,kp,ki)
    if size(m.Ared,1)==size(Q,1)
      Quse=Q
    else
      gv=N.gauge_vector(m);Quse=nullspace(reshape(gv/norm(gv),1,:))
    end
    A=Quse'*m.Ared*Quse
    F=eigen(A);lam=F.values;V=F.vectors;alpha=maximum(real.(lam))
    slow=sortperm(real.(lam);rev=true)[1:6]
    ws=unique(vcat(0.0,abs.(imag.(lam[slow])),extra_frequencies))
    betas=[minimum(svdvals(im*w*I-A-0.05I)) for w in ws]
    jb=argmin(betas);beta=betas[jb]
    t=collect(0.0:dt:horizon);T=0.5
    # Add both sides of the exact moving-window discontinuities.
    append!(t,[T-1e-9,T+1e-9,2T-1e-9,2T+1e-9]);sort!(unique!(t))
    nt=length(t);n=length(lam)
    function factors(shift)
      a=max.(t.-shift,0.0)
      H=zeros(ComplexF64,n,nt)
      for j in 1:n,k in 1:nt
        x=lam[j]*a[k]
        H[j,k]=abs(x)<1e-5 ? a[k]^2/2+lam[j]*a[k]^3/6+lam[j]^2*a[k]^4/24 :
          (expm1(x)-x)/lam[j]^2
      end
      H
    end
    H0=factors(0.0);H1=factors(T);H2=factors(2T)
    KF=(H0-H1)/T;KR=(H0-2H1+H2)/T^2
    DF=min.(t,T)/T;DR=(t.-2max.(t.-T,0.0).+max.(t.-2T,0.0))/T^2
    JF=Float64.(t.<T)/(2pi*T)
    JR=(Float64.(t.<T).-Float64.((t.>=T).&(t.<2T)))/(2pi*T^2)
    DV=-(m.Gy\LOADS);Bfull=m.B*DV;B=Quse'*Bfull
    C=-(ANG*(m.Gy\(m.C*m.Ared))*Quse)/(2pi)
    D=-(ANG*(m.Gy\(m.C*Bfull)))/(2pi)
    jump=ANG*DV;Z=V\B;CV=C*V
    peaks=zeros(3);rocofs=zeros(3);steady=zeros(3);pt=zeros(3);pb=zeros(Int,3)
    for b in 1:3
       residues=CV.*reshape(Z[:,b],1,:)
       freq=100 .* (real.(residues*KF)+D[:,b]*DF'+jump[:,b]*JF')
       roc=100 .* (real.(residues*KR)+D[:,b]*DR'+jump[:,b]*JR')
       ss=100abs.(D[:,b]-real.(CV*(Z[:,b]./lam)));steady[b]=maximum(ss)
       idx=argmax(abs.(freq));peaks[b]=max(steady[b],abs(freq[idx]));rocofs[b]=maximum(abs,roc)
       pb[b]=idx[1]+29;pt[b]=t[idx[2]]
       if steady[b]>abs(freq[idx]);pb[b]=argmax(ss)+29;pt[b]=Inf;end
    end
    # Positive g means violation. A 2% guard is imposed on the sampled radius.
    g=vcat((alpha+0.05)/0.05,max(-5.0,1-beta/(1.02CoreDesign.BETA_REQ)),
      peaks./0.5 .-1,rocofs./0.5 .-1)
    (;g,alpha,beta,omega=ws[jb],peaks,rocofs,steady,pt,pb,ep,kp,ki,A,lam,
      J=dot(CTX.power,ep))
end

if abspath(PROGRAM_FILE)==@__FILE__
  y=CoreDesign.encode(CTX,ORIGINAL["epsilon"],ORIGINAL["Kp"],ORIGINAL["Ki"])
  t=@elapsed v=oracle(y)
  println("ORACLE_CHECK J=",v.J," alpha=",v.alpha," F=",v.peaks," R=",v.rocofs," seconds=",t);flush(stdout)
  server=listen(ip"127.0.0.1",0);port=getsockname(server)[2]
  println("READY ",port);flush(stdout)
  sock=accept(server)
  while isopen(sock)
    line=readline(sock);line=="QUIT" && break
    try
      args=split(line,';');cmd=args[1];x=parse.(Float64,split(args[2],','))
      if cmd=="DYNJAC"
        out=vec(dynamic_jacobian(x))
      else
        vv=oracle(x;dt=cmd=="FINE" ? .01 : .05,horizon=cmd=="FINE" ? 90.0 : 40.0)
        out=vcat(vv.J,vv.alpha,vv.beta,vv.omega,vv.peaks,vv.rocofs,vv.steady,vv.g)
      end
      println(sock,join(out,','));flush(sock)
    catch ex
      println(sock,"ERR ",replace(sprint(showerror,ex),'\n'=>' '));flush(sock)
    end
  end
  close(sock);close(server)
end
