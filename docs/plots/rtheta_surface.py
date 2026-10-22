import matplotlib.pyplot as plt
from mpmath import j, rtheta, splot

tau = [[j, 0.5], [0.5, j]]
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surface = lambda x, y: abs(rtheta([x, y], tau))
splot(surface, [-1, 1], [-1, 1], points=35, keep_aspect=False, axes=ax, plot3d_kwargs={"cmap": "viridis"})
ax.set_zlabel(r"$|\theta(z\mid\tau)|$")
fig.savefig("rtheta_surface.png")
