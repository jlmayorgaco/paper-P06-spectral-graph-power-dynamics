module AlgebraicBoundary

using LinearAlgebra

export bilinear_gain_candidates, bilinear_resultant_coefficients

"Resultant polynomial coefficients in Kp from two bilinear real equations."
function bilinear_resultant_coefficients(c00,c10,c01,c11)
    a00,a10,a01,a11=real(c00),real(c10),real(c01),real(c11)
    b00,b10,b01,b11=imag(c00),imag(c10),imag(c01),imag(c11)
    # (b00+b10*x)(a01+a11*x)-(b01+b11*x)(a00+a10*x)
    return (c0=b00*a01-b01*a00,
            c1=b00*a11+b10*a01-b01*a10-b11*a00,
            c2=b10*a11-b11*a10)
end

"""
Solve Re(F)=Im(F)=0 for a demonstrated bilinear complex polynomial.
Returns all real gain pairs; bounds only classify candidates after solving.
"""
function bilinear_gain_candidates(c00,c10,c01,c11;
        kp_bounds=(-Inf,Inf),ki_bounds=(-Inf,Inf),tol=1e-10)
    q=bilinear_resultant_coefficients(c00,c10,c01,c11)
    scale=max(abs(q.c0),abs(q.c1),abs(q.c2),1.0)
    roots=Float64[]
    if abs(q.c2)<=tol*scale
        abs(q.c1)>tol*scale && push!(roots,-q.c0/q.c1)
    else
        disc=q.c1^2-4q.c2*q.c0
        disc>=-tol*scale^2 && begin
            d=sqrt(max(disc,0.0)); push!(roots,(-q.c1+d)/(2q.c2))
            d>tol*scale && push!(roots,(-q.c1-d)/(2q.c2))
        end
    end
    out=NamedTuple[]
    for kp in unique(roots)
        den=real(c01+c11*kp); num=real(c00+c10*kp)
        if abs(den)>tol*scale
            ki=-num/den
            F=c00+c10*kp+c01*ki+c11*kp*ki
            push!(out,(Kp=kp,Ki=ki,residual=abs(F),physical=kp_bounds[1]<=kp<=kp_bounds[2] &&
                ki_bounds[1]<=ki<=ki_bounds[2],exceptional=false))
        else
            # If the real equation loses Ki at this root, use the imaginary
            # equation and retain the candidate only if both equations agree.
            deni=imag(c01+c11*kp); numi=imag(c00+c10*kp)
            if abs(deni)>tol*scale
                ki=-numi/deni
                F=c00+c10*kp+c01*ki+c11*kp*ki
                push!(out,(Kp=kp,Ki=ki,residual=abs(F),physical=kp_bounds[1]<=kp<=kp_bounds[2] &&
                    ki_bounds[1]<=ki<=ki_bounds[2],exceptional=true))
            end
        end
    end
    return out
end

end
