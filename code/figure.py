"""Create the distance-power-one illustration using the proved limiting cost."""
from itertools import combinations
from math import sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"pdf.fonttype":42})
purple,green,ink = "#773F79","#15786A","#202332"
fig, (graph, ax) = plt.subplots(1, 2, figsize=(7.05,2.35),
                                gridspec_kw={"width_ratios":[1.05,1]})
fig.subplots_adjust(left=.025,right=.985,wspace=.32,bottom=.24,top=.84)
cmap=LinearSegmentedColormap.from_list("harmonic",[purple,"#AFA388",green])
path=np.c_[np.linspace(0,1,7), np.zeros(7)]
for offset, direction, colour in [(0,-1,purple),(1,1,green)]:
    angles=np.linspace(0,2*np.pi,8,endpoint=False)
    xy=np.c_[offset+direction*.30*(1-np.cos(angles)),.30*np.sin(angles)]
    for i,j in combinations(range(len(xy)),2):
        graph.plot(xy[[i,j],0],xy[[i,j],1],color=colour,alpha=.31,lw=.7,zorder=1)
    graph.scatter(xy[:,0],xy[:,1],s=24,color=colour,zorder=3)
for i in range(len(path)-1):
    graph.plot(path[i:i+2,0],path[i:i+2,1],lw=1.3,color=ink,zorder=1)
graph.scatter(path[:,0],path[:,1],c=cmap(np.linspace(0,1,len(path))),s=26,zorder=4)
graph.text(-.28,-.47,r"$y=+1$",ha="center",color=purple)
graph.text(1.28,-.47,r"$y=-1$",ha="center",color=green)
graph.text(.5,.20,r"$r$ edges",ha="center",color=ink)
graph.set(xlim=(-.7,1.7),ylim=(-.62,.55),aspect="equal")
graph.axis("off")
graph.set_title("An attaining construction",loc="left",fontsize=10,color=ink)
s=np.linspace(.6,3.5,260)
ax.plot(s,2/s,ls="--",c=purple,lw=1.2)
ax.plot(s,2*s/3,ls=":",c=green,lw=1.8)
ax.plot(s,2/s+2*s/3,c=ink,lw=1.8)
ax.plot([sqrt(3)],[4/sqrt(3)],"o",c=ink,ms=4)
ax.annotate(r"$4/\sqrt{3}$",xy=(sqrt(3),4/sqrt(3)),xytext=(1.92,2.85),
            fontsize=9,arrowprops={"arrowstyle":"-","color":ink,"lw":.7})
ax.text(.93,1.80,r"$2/s$",color=purple)
ax.text(2.80,1.57,r"$2s/3$",color=green)
ax.set_title("Scaled first-order cost",loc="left",fontsize=10,color=ink)
ax.set(xlim=(.6,3.5),ylim=(0,4.15),xlabel=r"Joining length $s=r/\sqrt{n}$",
       xticks=(1,sqrt(3),3),xticklabels=("1",r"$\sqrt{3}$","3"), yticks=(0,1,2,3,4))
ax.spines[["right","top"]].set_visible(False)
ax.spines[["bottom","left"]].set_color("#88838A")
ax.tick_params(axis="both",length=3,color="#88838A",labelsize=8)
fig.savefig(Path(__file__).resolve().parent.parent/"figures"/"distance-balance.pdf",bbox_inches="tight",pad_inches=.035)
