from close_design import *

def main():
    a=Model(ANCHOR);d=action(PAIR-ANCHOR);rows=[]
    for n in [128,256]:
        linear=0j;quadratic=0j
        for k in range(n):
            rr=.09*np.exp(2j*np.pi*(k+.5)/n);s=TARGET+rr;H=d[:,None]*a.T(s)
            linear-=np.trace(H)*rr/n;quadratic+=.5*np.trace(H@H)*rr/n
        zi=local(0,.02)[1];zj=local(7,.02)[1];joint=complex(*SELECTED['root']);nodal=zi+zj-TARGET
        for label,pred in [('first_order',TARGET+linear),('ordinary_quadratic',TARGET+linear+quadratic),('exact_singleton_additive',nodal)]:
            rows.append(dict(method=label,n=n,pred_real=pred.real,pred_imag=pred.imag,
                real_error=pred.real-joint.real,complex_error=abs(pred-joint),predicts_margin_pass=bool(pred.real<=-.05)))
        # C2 from same common box at the finer frozen numerical quadrature.
        dd=json.loads((OUT/'CONTOUR_DIAGNOSTICS.json').read_text())[1]
        pred=nodal+complex(*dd['C2'])
        rows.append(dict(method='singleton_plus_two_step',n=n,pred_real=pred.real,pred_imag=pred.imag,
            real_error=pred.real-joint.real,complex_error=abs(pred-joint),predicts_margin_pass=bool(pred.real<=-.05)))
    pd.DataFrame(rows).to_csv(OUT/'TABLE_09_PREDICTOR_COMPARISON.csv',index=False)
    # Independent centered differences validate directional coefficients only.
    hh=1e-3;zp=Model(ANCHOR+hh*(PAIR-ANCHOR)).refine(TARGET)[0];zm=Model(ANCHOR-hh*(PAIR-ANCHOR)).refine(TARGET)[0]
    save('PREDICTOR_DERIVATIVE_CHECK.json',{'linear_error':abs((zp-zm)/(2*hh)-linear),
        'quadratic_coefficient_error':abs((zp+zm-2*TARGET)/(2*hh*hh)-quadratic),
        'series_finite_remainder':'NOT_VALID_ON_SELECTED_COMMON_CONTOUR','ranked_methods':'counterfactual surrogates, not reproductions of published algorithms'})
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
