# Mathematical one-sided limits only. No zero-rated device is accepted as an
# actual architecture. These limits test whether nearby insertions are feasible.
function adjacent_exact_limits(d)
    old=L.N.descriptor(L.CTX,d["rho"],d["Kp"],d["Ki"])
    go=L.N.gauge_vector(old);Qo=nullspace(reshape(go/norm(go),1,:));Ao=Qo'*old.Ared*Qo
    rows=NamedTuple[]
    for bus in setdiff(30:39,Int.(d["support"]))
        rho=copy(d["rho"]);rho[bus-29]=.999
        m=L.N.descriptor(L.CTX,rho,d["Kp"],d["Ki"])
        C=copy(m.C);D=copy(m.D);yi=2bus-1:2bus
        for (block,row) in zip(m.blocks,eachrow(m.state_map))
            block.bus==bus || continue
            xi=Int(row.first):Int(row.last)
            share=block.kind=="SG" ? 0. : 1.
            C[yi,xi].=share.*block.J.C
            D[yi,yi].+=(share-block.share).*block.J.D
        end
        Ar=m.A-m.B*((L.CTX.net.y_static+D)\C)
        gn=L.N.gauge_vector(m);Qn=nullspace(reshape(gn/norm(gn),1,:));An=Qn'*Ar*Qn
        S=zeros(size(old.A,1),size(m.A,1))
        for row in eachrow(old.state_map)
            nr=only(eachrow(m.state_map[(m.state_map.bus.==row.bus).&(m.state_map.kind.==row.kind),:]))
            S[Int(row.first):Int(row.last),Int(nr.first):Int(nr.last)].=I(Int(row.last-row.first+1))
        end
        P=Qo'*S*Qn
        # A pointwise inverse norm lower bound provides a beta UPPER bound,
        # sufficient to prove violation. Power iteration only needs a witness.
        Rb=inv(-BigFloat.(An+.05I));y=BigFloat.(ones(size(An,1)));y/=norm(y)
        for _ in 1:80;y=Rb'*(Rb*y);y/=norm(y);end
        beta_upper=Float64(1/norm(Rb*y))
        push!(rows,(;bus,alpha_limit=maximum(real.(eigvals(An))),beta_zero_upper=beta_upper,
            beta_requirement=L.BETA,robust_violation=L.BETA-beta_upper,
            insertion_limit_infeasible=beta_upper<L.BETA,
            coisometry_residual=norm(P*P'-I),intertwining_residual=norm(P*An-Ao*P),
            gauge_residual=norm(Ar*gn)/norm(Ar)/norm(gn),
            inverse_residual_big=Float64(norm(I+BigFloat.(An+.05I)*Rb))))
    end
    CSV.write(joinpath(L.OUT,"ADJACENT_EXACT_LIMITS.csv"),DataFrame(rows))
    println(DataFrame(rows));flush(stdout)
    rows
end
