using CSV, DataFrames, Graphs, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__,"..","..","src","bnd_model_audit_m","PDCasesM.jl"))
using .PDCasesM

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_M","tables")

function evaluate(label,nw)
    eq=PDCasesM.PD39.initialize_equilibrium(nw;sparse=false,check=:error)
    λ=jacobian_eigenvals(eq.state)
    gauges=sort(abs.(λ))[1:2]
    finite=λ[abs.(λ).>1e-8]
    α=maximum(real.(finite))
    v30=Float64(eq.state[VIndex(30,:busbar₊Vbase)])
    v36=Float64(eq.state[VIndex(36,:busbar₊Vbase)])
    return (;label,alpha_excluding_numeric_zeros=α,raw_alpha=maximum(real.(λ)),
        Vbase30_kV=v30,Vbase36_kV=v36,smallest_pole_magnitudes=string(gauges),
        equilibrium_residual=PDCasesM.PD39._eq_residual(eq.state),spectrum=λ)
end

function main()
    println("VBASE_ABLATION_BUILD");flush(stdout)
    original=PDCasesM.build_case("ExpK_nominal")
    vertices,edges=PDCasesM.PD39.PD39Model.copy_network_components(original)
    data=PDCasesM.PD39.PD39Model.ieee39_data()
    for bus in 1:39
        r=data.bus[findfirst(==(bus),data.bus.bus),:]
        set_default!(vertices[bus],:busbar₊Vbase,Float64(r.base_kv))
    end
    corrected=Network(vertices,edges)
    set_jac_prototype!(corrected)
    a=evaluate("original_ExpK_builder",original)
    b=evaluate("CSV_Vbase_corrected",corrected)
    λa=sort(ComplexF64.(a.spectrum);by=x->(real(x),imag(x)))
    λb=sort(ComplexF64.(b.spectrum);by=x->(real(x),imag(x)))
    rows=[(;label=a.label,alpha_excluding_numeric_zeros=a.alpha_excluding_numeric_zeros,
        raw_alpha=a.raw_alpha,Vbase30_kV=a.Vbase30_kV,Vbase36_kV=a.Vbase36_kV,
        equilibrium_residual=a.equilibrium_residual,max_spectrum_difference=maximum(abs.(λa.-λb))),
      (;label=b.label,alpha_excluding_numeric_zeros=b.alpha_excluding_numeric_zeros,
        raw_alpha=b.raw_alpha,Vbase30_kV=b.Vbase30_kV,Vbase36_kV=b.Vbase36_kV,
        equilibrium_residual=b.equilibrium_residual,max_spectrum_difference=maximum(abs.(λa.-λb)))]
    CSV.write(joinpath(OUT,"TABLE_M13_voltage_base_ablation.csv"),DataFrame(rows))
    println("VBASE_ABLATION_DONE alpha_original=",a.alpha_excluding_numeric_zeros,
        " alpha_corrected=",b.alpha_excluding_numeric_zeros,
        " max_delta=",maximum(abs.(λa.-λb)));flush(stdout)
end

main()
