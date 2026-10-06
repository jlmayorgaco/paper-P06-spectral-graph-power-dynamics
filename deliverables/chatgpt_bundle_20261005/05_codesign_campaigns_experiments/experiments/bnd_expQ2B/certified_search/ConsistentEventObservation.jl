using ForwardDiff
# At a parameter jump, interpolation may combine pre-jump algebraic voltages
# with post-jump parameters. Reconstruct ONLY the algebraic right limit while
# holding all differential states fixed. No integration, control or trim change.
function consistent_event_observation(sol,t)
    s=NWState(sol,t);nw=NetworkDynamics.extract_nw(s)
    u=copy(uflat(s));p=copy(pflat(s));mass=nw.mass_matrix
    md=mass isa UniformScaling ? ones(length(u)) : diag(mass)
    ai=findall(iszero,md);di=findall(!iszero,md);ud=copy(u[di])
    function algebraic(y)
        xx=eltype(y).(u);xx[ai]=y;f=similar(xx);nw(f,xx,p,t);f[ai]
    end
    for _ in 1:8
        y=u[ai];f=algebraic(y)
        norm(f,Inf)<1e-10 && break
        J=ForwardDiff.jacobian(algebraic,y);u[ai]-=J\f
    end
    err=norm(algebraic(u[ai]),Inf)
    err<1e-9 || error("Post-event algebraic observation residual $err")
    u[di]==ud || error("Differential state changed during observation")
    NWState(nw,u,p,t)
end
