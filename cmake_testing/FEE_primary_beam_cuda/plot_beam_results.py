import matplotlib.pyplot as plt
import numpy as np


def make_ax_do_plot(fig, az, za, plot_ind, data, numx=2, numy=4,
                    s=15, vmax=False, alpha=0.5 ):
    ax = fig.add_subplot(numy,numx,plot_ind,projection='polar')

    if vmax:
        ax.scatter(az, za, c=data, s=s, alpha=alpha, vmax=vmax)
    else:
        ax.scatter(az, za, c=data, s=s, alpha=alpha)
    ax.set_ylim(0,np.pi/2)
    ax.set_yticklabels([])

az, za, g1r, g1i, D1r, D1i, D2r, D2i, g2r, g2i = np.loadtxt('beam_values_out.txt', unpack=True)


# fig = plt.figure(figsize=(6,12))
#
# data_sets = [g1r, g1i, D1r, D1i, D2r, D2i, g2r, g2i]
#
# for ind, data in enumerate(data_sets):
#     make_ax_do_plot(fig, az, za, ind + 1, data)
#
# plt.tight_layout()
# fig.savefig('beam_values.png',bbox_inches='tight')
# plt.close()


g1 = g1r + 1j*g1i
D1 = D1r + 1j*D1i
D2 = D2r + 1j*D2i
g2 = g2r + 1j*g2i

XX = g1*np.conjugate(g1) + D1*np.conjugate(D1)
YY = g2*np.conjugate(g2) + D2*np.conjugate(D2)



fig = plt.figure(figsize=(10,10))

data_sets = [np.real(XX), np.imag(XX), np.real(YY), np.imag(YY)]

for ind, data in enumerate(data_sets):
    make_ax_do_plot(fig, az, za, ind + 1, data, numx=2, numy=2,
                    s=30, vmax=0.3, alpha=0.7)

plt.tight_layout()
fig.savefig('beam_values_XX-YY.png',bbox_inches='tight')
plt.close()
