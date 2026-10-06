module PortReduction

using LinearAlgebra
using ..AnalyticDeviceModel: port_transfer

export device_admittance, rated_parallel_mix, affine_rho_residual

device_admittance(A,B,C,D,s) = port_transfer(A,B,C,D,s)

"""
Parallel-current aggregation for device port matrices already referred to the
common system base. The caller must perform each device's MVA-base conversion
before calling this helper; it does not declare an SG or GFL model physically
valid by itself.
"""
function rated_parallel_mix(Ysg::AbstractMatrix,Ygfl::AbstractMatrix,rho::Real)
    0<=rho<=1 || throw(ArgumentError("rho must be in [0,1]"))
    size(Ysg)==size(Ygfl) || throw(DimensionMismatch("device port maps must share a port convention"))
    return (1-rho).*Ysg .+ rho.*Ygfl
end

function affine_rho_residual(Ysg,Ygfl,rho)
    Ymix=rated_parallel_mix(Ysg,Ygfl,rho)
    return norm(Ymix-(Ysg+rho.*(Ygfl-Ysg)))/max(norm(Ymix),eps(Float64))
end

end
