using Sockets, LinearAlgebra, SHA, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const CANDIDATE=joinpath(ROOT,"reports","experiment_N","Z_N_NOMINAL_FINAL.toml")
bytes2hex(sha256(read(CANDIDATE)))==strip(read(CANDIDATE*".sha256",String)) ||
    error("generic benchmark requires frozen candidate")
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
using .PDExactDesignN
ctx=design_context(ROOT)
server=listen(ip"127.0.0.1",32761)
println("FALSIFICATION_SERVER_READY");flush(stdout)
sock=accept(server)
count=0
while isopen(sock)
    line=try readline(sock) catch; break end
    line=="QUIT" && break
    try
        z=parse.(Float64,split(line,','))
        length(z)==21 || error("21 variables required")
        e=z[1];kp=z[2:11];ki=z[12:21]
        0<e<=1 || error("eps outside fixed-support interval")
        rho=ones(10);rho[9]=1-e
        sp=spectrum(ctx,rho,kp,ki)
        d=simple_mode_sensitivities(ctx,rho,kp,ki)
        values=vcat(sp.alpha,-real(d.rho[9]),real.(d.Kp),real.(d.Ki))
        println(sock,join(values,','));flush(sock)
        global count=count+1
    catch ex
        println(sock,"ERR:",replace(sprint(showerror,ex),'\n'=>' '));flush(sock)
    end
end
close(sock);close(server)
println("FALSIFICATION_SERVER_DONE evaluations=",count);flush(stdout)
