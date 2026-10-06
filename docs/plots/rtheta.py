import matplotlib.pyplot as plt
from mpmath import j, plot, re, rtheta

tau = [[j, -0.5], [-0.5, j]]
curves = [
    lambda x: re(rtheta([x, x/2], tau)),
    lambda x: re(rtheta([x, 2*x], tau)),
]
fig, ax = plt.subplots()
plot(curves, [-2, 2], axes=ax)
ax.legend([r"$z=(x,x/2)$", r"$z=(x,2x)$"])
fig.savefig("rtheta.png")
