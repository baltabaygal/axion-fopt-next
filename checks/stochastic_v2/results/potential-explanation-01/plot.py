"""Plot the exact bare potential used in the first bubble pilot."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
s=np.linspace(-.2,1.2,1000)
u=lambda x:.5*x*x*(1-x)**2+.06*(1-4*x**3+3*x**4)
b=2/(4+24*.06)
fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
ax.plot(s,u(s),color='#315c9b',lw=2.5)
for x,label,text in [(0,'False vacuum\nS = 0, U = 0.06',(.03,.102)),(b,'Barrier',(.52,.094)),(1,'True vacuum\nS = 1, U = 0',(.9,.028))]:
    ax.scatter([x],[u(x)],color='#315c9b',s=35,zorder=3)
    ax.annotate(label,(x,u(x)),xytext=text,ha='center',arrowprops=dict(arrowstyle='-',color='gray'))
ax.axhline(0,color='gray',lw=.6)
ax.set(xlabel='Field value S (not position)',ylabel='Potential energy density U(S)',xlim=(-.2,1.2),ylim=(-.008,.12),title='The transition-field potential used in our first experiment')
ax.spines[['top','right']].set_visible(False)
fig.savefig(Path(__file__).with_name('potential.png'),dpi=170)
