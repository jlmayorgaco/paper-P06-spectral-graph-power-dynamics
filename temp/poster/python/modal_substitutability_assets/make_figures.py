"""
POSTER FIGURE GENERATOR
Figure 1 (dominant): 6-bus graph, nodes colored by modal substitutability s_i.
  green = topology-cheap (high s), red = inertia-inevitable (low s).
Figure 2: paired bar chart - inertia shift (down) vs line shift (up) on nu_c.
"""
import numpy as np, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

edges=[(0,1,1.2),(0,2,0.8),(1,2,1.5),(1,3,0.6),(2,4,1.0),(3,4,0.9),(3,5,1.3),(4,5,0.7)]
N=6
def build_L(e,N):
    L=np.zeros((N,N))
    for(i,j,w) in e:L[i,j]-=w;L[j,i]-=w;L[i,i]+=w;L[j,j]+=w
    return L
L=build_L(edges,N);M=np.diag([2.0,0.5,3.0,0.8,1.5,0.4])
ms=np.diag(1/np.sqrt(np.diag(M)));Lt0=ms@L@ms
nu,Q=np.linalg.eigh(Lt0);crit=np.argmax(nu);q=Q[:,crit]
subst=np.array([1-abs((0.5*(L@np.eye(N)[i]))@q)/(np.linalg.norm(0.5*(L@np.eye(N)[i]))+1e-12) for i in range(N)])

# node layout (hand-placed for the 6-bus)
pos={0:(0,1),1:(1,1.6),2:(1,0.4),3:(2,1.6),4:(2,0.4),5:(3,1)}

# FIGURE 1 - substitutability graph
fig,ax=plt.subplots(figsize=(8,5))
for (i,j,w) in edges:
    x=[pos[i][0],pos[j][0]];y=[pos[i][1],pos[j][1]]
    ax.plot(x,y,color="#888",lw=1+2*w,zorder=1,alpha=0.6)
cmap=plt.cm.RdYlGn
for n in range(N):
    c=cmap(subst[n])
    ax.scatter(*pos[n],s=2200,color=c,edgecolors="k",lw=2,zorder=3)
    ax.text(*pos[n],f"{n}\ns={subst[n]:.2f}",ha="center",va="center",fontsize=11,fontweight="bold",zorder=4)
ax.set_title("Modal substitutability of virtual inertia\n(green = replace with lines · red = inertia/SynCon needed)",fontsize=13,fontweight="bold")
sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(0,1));sm.set_array([])
cb=plt.colorbar(sm,ax=ax,fraction=0.04);cb.set_label("substitutable by topology  →",fontsize=11)
ax.axis("off");plt.tight_layout();plt.savefig("/home/claude/fig1_substitutability.png",dpi=200,bbox_inches="tight")
plt.savefig("/home/claude/fig1_substitutability.svg",bbox_inches="tight");plt.close()

# FIGURE 2 - opposite-direction levers
fig,ax=plt.subplots(figsize=(8,4.5))
inertia=[]
for i in range(N):
    Ei=np.zeros((N,N));Ei[i,i]=1;dLt=-0.5*(Ei@Lt0+Lt0@Ei);inertia.append(q@dLt@q)
lines=[]
lbls=[]
for (i,j,w) in edges:
    b=np.zeros(N);b[i]=1;b[j]=-1;dLt=ms@np.outer(b,b)@ms;lines.append(q@dLt@q);lbls.append(f"{i}-{j}")
x1=np.arange(N)
ax.bar(x1,inertia,color="#c0392b",label="inertia at node (lowers $\\nu_c$)")
ax.set_xticks(x1);ax.set_xticklabels([f"n{i}" for i in range(N)])
ax2=ax.twiny();x2=np.arange(len(lines))+0.0
ax.axhline(0,color="k",lw=0.8)
# overlay line bars on a second axis region (just show as separate group)
ax.bar(np.arange(N,N+len(lines)),lines,color="#27ae60",label="line reinforcement (raises $\\nu_c$)")
ax.set_xticks(list(range(N))+list(range(N,N+len(lines))))
ax.set_xticklabels([f"n{i}" for i in range(N)]+lbls,rotation=45,fontsize=8)
ax2.axis("off")
ax.set_ylabel("first-order shift of $\\nu_c$",fontsize=11)
ax.set_title("The two levers move the critical mode in OPPOSITE directions",fontsize=12,fontweight="bold")
ax.legend(fontsize=10);plt.tight_layout()
plt.savefig("/home/claude/fig2_opposite_levers.png",dpi=200,bbox_inches="tight")
plt.savefig("/home/claude/fig2_opposite_levers.svg",bbox_inches="tight");plt.close()
print("figures written: fig1_substitutability, fig2_opposite_levers (png+svg)")
print("substitutability:",np.round(subst,3))
