"""Extract the frozen explicit Julia model functions without rerunning exports."""
from pathlib import Path
D=Path(__file__).resolve().parent
s=(D/"export_coupled_dq_model.jl").read_text(encoding="utf-8")
head=s[:s.index("\nz=zeros(nx)")]
(D/"CoupledDQReference.jl").write_text('module CoupledDQReference\n'+head+'\nend\n',encoding="utf-8")
print("Prepared Julia reference with explicit full_rhs/global_rhs; no report export side effects.")
