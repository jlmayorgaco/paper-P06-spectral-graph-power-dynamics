import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from t8_taylor_and_second_design import *
for rho in (0.90, 0.80):
    p = P0.copy(); p[:10] = rho
    sel = select_node(p)[0]
    print(rho, repr(sel[0].real), repr(sel[0].imag), sel[0].imag/6.283185307, flush=True)
