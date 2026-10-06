using LinearAlgebra

"Build a noncommutative graph-feature basis from the frozen port operator."
function graph_features(net,buses,dispatch,gamma)
    n=length(buses)
    P=zeros(Float64,size(net.y_static,1),2n)
    for (j,b) in enumerate(buses); P[2b-1:2b,2j-1:2j].=I(2); end
    G=transpose(P)*(net.y_static\P)
    Y=inv(G)
    LG=(Y+transpose(Y))/2
    LB=(Y-transpose(Y))/2
    x=repeat(Float64.(dispatch)./max(norm(dispatch),eps()),inner=2)
    words=[x,LG*x,LB*x,LG^2*x,LB^2*x,LG*LB*x,LB*LG*x]
    raw=zeros(Float64,n,length(words))
    for (k,w) in enumerate(words), i in 1:n
        raw[i,k]=norm(w[2i-1:2i])
    end
    strengths=[norm(Y[2i-1:2i,:]) for i in 1:n]
    auth=repeat(Float64.(gamma)./max(norm(gamma),eps()),inner=1)
    raw=hcat(raw,dispatch./max(norm(dispatch),eps()),strengths./max(norm(strengths),eps()),auth)
    F=qr(raw)
    Q=Matrix(F.Q)[:,1:min(size(raw)...)]
    rankQ=count(>(1e-10),diag(F.R)[1:min(size(F.R)...)]); Q=Q[:,1:rankQ]
    (;LG,LB,raw_features=raw,basis=Q,rank=rankQ,
      names=["I*x","LG*x","LB*x","LG2*x","LB2*x","LG_LB*x","LB_LG*x",
             "dispatch","port_strength","SG_authority"],
      simultaneous_diagonalization_claim=false)
end
