import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.lines as mlines
import matplotlib.font_manager
import seaborn as sns
import pandas as pd
M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RefL0100N1504.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))
#M_h, m_s, t_i, t_q, sfr = np.loadtxt("data_plot_RecalL0025N0752.txt", unpack=True, usecols=( 2, 3, 4, 5, 6))


logmassmin = 9
logmassmax = 12
dex = 0.25
bin_num = (logmassmax - logmassmin)/dex + 1
bin_range = np.linspace(9, 12, bin_num)

keep_these = np.where(t_q != 14.134423953091941) #
and_these = np.where((t_i<=t_q))
q_before_i = np.where((t_i>=t_q))
unquenched = np.where(t_q == 14.134423953091941)
print(type(and_these))
# print(len(keep_these))
#print(and_these[0])
nuq = len(m_s)- len(m_s[keep_these]) + len(m_s[and_these])

#print(m_s[keep_these])
a =set((and_these[0]))
b = set((keep_these[0]))
c = b.intersection(a)
# c = np.array(list(c))


m9 = set((np.where((m_s >= 10**9) & (m_s <10**9.5))[0]))
m95 = set((np.where((m_s >= 10**9.5) & (m_s <10**10))[0]))
m10 = set((np.where((m_s >= 10**10) & (m_s <10**10.5))[0]))
m105 = set((np.where((m_s >= 10**10.5) & (m_s <10**11))[0]))
m11 = set((np.where((m_s >= 10**11))[0]))


c9 = c.intersection(m9)
c95 = c.intersection(m95)
c10= c.intersection(m10)
c105= c.intersection(m105)
c11= c.intersection(m11)

c9 = np.array(list(c9))
c95 = np.array(list(c95))
c10 = np.array(list(c10))
c105 = np.array(list(c105))
c11 = np.array(list(c11))

c = np.array(list(c))
import matplotlib
matplotlib.matplotlib_fname()

#font = {'family' : 'monospace',
#        'size'   : '12'}

font = {'size'   : '8'}
plt.rc('text', usetex=True)
print("quenched after infall: {0} \nquenched before infall: {1}\n unquenched: {2}\n total {3}".format(len(c),
                                                                                len(m_s[q_before_i]),len(m_s[unquenched]),len(m_s)))

quench_timescale = t_q[c] - t_i[c]
# bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=[np.log10(1e9),
#                             np.log10(10**9.5) ,np.log10(1e10),np.log10(10**10.5) ,np.log10(1e11), np.log10(10**11.5) ,np.log10(1e12)])

bin_means, bin_edges, binnumber = stats.binned_statistic( np.log10(m_s[c]),quench_timescale, statistic='median', bins=bin_range)

print(bin_means,bin_edges, binnumber)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2

#plt.rc('text', usetex=False)

# x = np.log10(m_s[c])
# bins = np.array([np.log10(1e9),np.log(10**9.5) ,np.log(1e10),np.log(10**10.5) ,np.log(1e11), np.log(10**11.5) ,np.log(1e12)])
# digitized = np.digitize(x, bins)
# print(digitized)
# print(x)
# print(np.log10(10))
# bin_means = [x[digitized == i].mean() for i in range(1, len(bins))]
# print(bin_means)
#fig = plt.figure(figsize=[8,10])
#data = pd.DataFrame(np.array([t_q[c]-t_i[c], t_i[c]]), columns=['x','y'])

#data = np.random.multivariate_normal([0, 0], [[5, 2], [2, 2]], size=2000)
#print(data)
#data = zip(t_q[c]-t_i[c], t_i[c])
#data = pd.DataFrame(data, columns=['x', 'y'])





fig = plt.figure(figsize=(1.5*6, 1.5*8.56))
grid = plt.GridSpec(11*2, 7*2)
grid.update(wspace=0., hspace=0)
main_ax0 = fig.add_subplot(grid[0:2*2, 1*2:3*2])
y_hist0 = fig.add_subplot(grid[0:2*2, 0:2], xticklabels=[], sharey=main_ax0)
x_hist0 = fig.add_subplot(grid[2*2:2*3, 1*2:3*2], yticklabels=[], sharex=main_ax0)

main_ax1 = fig.add_subplot(grid[0:2*2, 1*2+8:3*2+8])
y_hist1 = fig.add_subplot(grid[0:2*2, 0+8:2+8], sharey=main_ax1)
x_hist1 = fig.add_subplot(grid[2*2:2*3, 1*2+8:3*2+8], sharex=main_ax1)
x_hist1.set_ylim(0,115)

main_ax2 = fig.add_subplot(grid[0+8-1:2*2+8-1, 1*2:3*2])
y_hist2 = fig.add_subplot(grid[0+8-1:2*2+8-1, 0:2], sharey=main_ax2)
x_hist2 = fig.add_subplot(grid[2*2+8-1:2*3+8-1, 1*2:3*2], sharex=main_ax2)

main_ax3 = fig.add_subplot(grid[0+8-1:2*2+8-1, 1*2+8:3*2+8])
y_hist3 = fig.add_subplot(grid[0+8-1:2*2+8-1, 0+8:2+8], sharey=main_ax3)
x_hist3 = fig.add_subplot(grid[2*2+8-1:2*3+8-1, 1*2+8:3*2+8], sharex=main_ax3)

main_ax4 = fig.add_subplot(grid[0+8*2-2:2*2+8*2-2, 1*2:3*2])
y_hist4 = fig.add_subplot(grid[0+8*2-2:2*2+8*2-2, 0:2], sharey=main_ax4)
x_hist4 = fig.add_subplot(grid[2*2+2*8-2:2*3+8*2-2, 1*2:3*2], sharex=main_ax4)

main_ax5 = fig.add_subplot(grid[0+8*2-2:2*2+8*2-2, 1*2+8:3*2+8])
y_hist5 = fig.add_subplot(grid[0+8*2-2:2*2+8*2-2, 0+8:2+8], sharey=main_ax5)
x_hist5 = fig.add_subplot(grid[2*2+8*2-2:2*3+8*2-2, 1*2+8:3*2+8], sharex=main_ax5)


main_ax = np.array([main_ax0,main_ax1,main_ax2, main_ax3, main_ax4, main_ax5])
y_hist = np.array([y_hist0,y_hist1,y_hist2,y_hist3,y_hist4,y_hist5])
x_hist = np.array([x_hist0,x_hist1,x_hist2,x_hist3,x_hist4,x_hist5])

major_xticks = np.linspace(0, 14, 8)

# Set major ticks for y axis
major_yticks = np.linspace(0, 14, 8)

# I want minor ticks for x axis
minor_xticks = np.linspace(0, 14, 15)

# I want minor ticks for y axis
minor_yticks = np.linspace(0, 14, 15)


tau_q_d=[t_q[c9]-t_i[c9],t_q[c95]-t_i[c95],t_q[c10]-t_i[c10],t_q[c105]-t_i[c105],t_q[c11]-t_i[c11]]
t_i_d = [t_i[c9],t_i[c95],t_i[c10],t_i[c105],t_i[c11]]
# Specify tick label size
label_m = [r'$M_{\star} \rm /M_{\odot} =[9,9.5]$',r'$ M_{\star} \rm /M_{\odot} =[9.5,10]$', r'$M_{\star} \rm /M_{\odot}=[10,10.5]$', r'$M_{\star} \rm /M_{\odot}=[10.5,11]$',r'$M_{\star} \rm /M_{\odot}=[11,12]$']
empty_ax = fig.add_subplot(grid[2*2:2*3, 0:2], xticklabels=[], sharey=y_hist0,sharex=x_hist0)
empty_ax.axis("off")
empty_ax.hist([[-10],[-10],[-10],[-10],[-10]],color=['#032F64','#034698','#006CBB','#28A7EA','#45BDEE'],label=label_m)
empty_ax.legend(frameon=False)
# Suppress minor tick labels

for i in range(6):

    main_ax[i].set_xticks(major_xticks)
    main_ax[i].set_xticks(minor_yticks, minor = True)

    main_ax[i].set_yticks(major_xticks)
    main_ax[i].set_yticks(minor_yticks, minor = True)
    main_ax[i].tick_params(labelbottom=False,labelleft=False)

    main_ax[i].set_ylim(0,14)
    main_ax[i].set_xlim((0,14))
    y_hist[i].set_ylabel(r"$ \rm \tau_q [Gyr]$")
    x_hist[i].set_xlabel(r"$\rm t_i \ [Gyr]$")


main_ax[0].scatter(t_i[c] ,t_q[c]-t_i[c], marker = 'o',color='#034698', s=9, alpha=0.6,label = r'$ M_{\star} \rm /M_{\odot} =[9,12]$')
leg = main_ax[0].legend()
leg.get_frame().set_edgecolor('k')
# histogram on the attached axes
x_hist[0].hist(t_i_d, 28, histtype='step',stacked=True,
            orientation='vertical', color=['#032F64','#034698','#006CBB','#28A7EA','#45BDEE'])
#x_hist.invert_yaxis()

y_hist[0].hist(tau_q_d, 28,stacked=True, histtype='step',
            orientation='horizontal',color=['#032F64','#034698','#006CBB','#28A7EA','#45BDEE'])
#y_hist.invert_xaxis()

plt.tight_layout()



#for j in range(1,6):


for j in range(5):
    main_ax[j+1].scatter(t_i_d[j], tau_q_d[j],marker='o', color='#034698', s=9, alpha=0.7,
                       label=label_m[j])

    x_hist[j+1].hist(t_i_d[j], 28, histtype='step', stacked=True,
                   orientation='vertical', color='#034698')
    # x_hist.invert_yaxis()

    y_hist[j+1].hist(tau_q_d[j], 28, stacked=True, histtype='step',
                   orientation='horizontal', color='#034698')
    leg = main_ax[j+1].legend()
    leg.get_frame().set_edgecolor('k')

plt.savefig('plots20.06.19/tau_q_vs_t_i_tight_layout.pdf', format='pdf', dpi=900, pad_inches=0.1, bbox_inches='tight')
plt.show()

