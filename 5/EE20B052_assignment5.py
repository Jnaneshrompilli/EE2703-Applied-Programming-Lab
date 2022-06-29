"""-------------------------------------------------------------------

Author: Jnaneswara Rao Rompilli (EE20B052)
Date: 17-2-22, 22:22
Description: Solving currents using Laplace Equation

---------------------------------------------------------------------"""

from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
from sys import argv

Nx = 25  # size laong x
Ny = 25  # size along y
radius = 8  # radius of central lead
Niter = 1500  # number of iterations to perform

argc = len(argv)

# Command line arguments
if argc == 2:
    Nx = argv[1]

elif argc == 3:
    Nx = argv[1]
    Ny = argv[2]

elif argc == 4:
    Nx = argv[1]
    Ny = argv[2]
    radius = argv[3]

elif argc == 5:
    Nx = argv[1]
    Ny = argv[2]
    radius = argv[3]
    Niter = argv[4]

# Function to fit exponential curve using lstsq
def fit_curve(Niter, start, errors):
    xfir = np.full((Niter - start, 1), fill_value=1)
    xsec = np.arange(start, Niter).reshape((Niter - start, 1))
    M = np.hstack((xfir, xsec))
    y2 = log(errors[start:])
    coeff2 = lstsq(M, y2, rcond=None)[0]
    return exp(coeff2[0]) * exp(coeff2[1] * xsec)


# Intialize required variables
phi = np.zeros((Ny, Nx))
x = np.linspace(-0.5, 0.5, Nx)
y = np.linspace(-0.5, 0.5, Ny)

Y, X = np.meshgrid(y, x)

# Assign potential = 1
ii = np.where(X * X + Y * Y <= (0.35 * 0.35))
phi[ii] = 1

# Plot contour of V
plt.figure(1, figsize=(7, 7))
plt.plot(x[ii[0]], y[ii[1]], "or", label="V = 1")
plt.title("Contour plot of potential", fontsize=16)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.xlabel("x", fontsize=14)
plt.ylabel("y", fontsize=14)
grid()
plt.legend()
# plt.savefig("1.png")

oldphi = np.zeros((Ny, Nx))
errors = np.zeros((Niter, 1))

# Updating the potential
for k in range(Niter):
    oldphi = phi.copy()
    new = oldphi[2:, 1:-1] + oldphi[1:-1, 0:-2] + oldphi[1:-1, 2:] + oldphi[0:-2, 1:-1]
    phi[1:-1, 1:-1] = 0.25 * new

    # Vectorizing the code
    phi[1:-1, 0] = phi[1:-1, 1]
    phi[1:-1, -1] = phi[1:-1, -2]
    phi[0, 1:-1] = phi[1, 1:-1]
    phi[ii] = 1.0

    errors[k] = (abs(phi - oldphi)).max()


# logy = logA + Bx
yfit1 = fit_curve(Niter, 0, errors)
yfit2 = fit_curve(Niter, 500, errors)

k = np.arange(0, Niter)

# Plotting errors in semilog scale
plt.figure(2, figsize=(7, 7))
plt.semilogy(k, errors, "go", markersize=2, label="errors")
plt.xlabel("Nth iteration", fontsize=14)
plt.ylabel("Error", fontsize=14)
plt.title("In phi matrix (logy scale)", fontsize=16)
grid()
plt.legend()
# plt.savefig("2.png")

# Plotting errors in loglog scale
plt.figure(3, figsize=(7, 7))
plt.loglog(k, errors, "go", markersize=2, label="errors")
plt.xlabel("Nth iteration", fontsize=14)
plt.ylabel("Error", fontsize=14)
plt.title("In phi matrix (loglog scale)", fontsize=16)
grid()
plt.legend()
# plt.savefig("3.png")

# Plotting errors after 500th iteration in loglog scale
plt.figure(4, figsize=(7, 7))
plt.semilogy(k[500:], errors[500:], "go", markersize=2, label="errors")
plt.xlabel("Nth iteration", fontsize=14)
plt.ylabel("Error", fontsize=14)
plt.title("In phi matrix (loglog scale)", fontsize=16)
grid()
plt.legend()
# plt.savefig("4.png")

# Plotting errors and fitted curves in loglog scale
figure(5, figsize=(7, 7))
plt.loglog(k, errors, "+r", markersize="2", label="errors")
plt.loglog(k, yfit1, linewidth=2, label="fit1")
plt.loglog(k[500:], yfit2, linewidth=2, label="fit2", color="green")
plt.title("Exponential fitting", fontsize=16)
plt.xlabel("Nth iteration", fontsize=14)
plt.ylabel("Error", fontsize=14)
grid()
plt.legend()
# plt.savefig("5.png")

# Plotting Potential
plt.figure(6, figsize=(7, 7))
plt.contourf(Y, X[::-1], phi)
plt.plot(x[ii[0]], y[ii[1]], "or", label="V = 1")
plt.title("Contour plot of V", fontsize=16)
plt.xlabel("X", fontsize=14)
plt.ylabel("Y", fontsize=14)
grid()
plt.legend()
# plt.savefig("6.png")

# Currents
Jx, Jy = (
    1 / 2 * (phi[1:-1, 0:-2] - phi[1:-1, 2:]),
    1 / 2 * (phi[:-2, 1:-1] - phi[2:, 1:-1]),
)

# Vector plot of currents
figure(7, figsize=(7, 7))
plt.quiver(Y[1:-1, 1:-1], -X[1:-1, 1:-1], -Jx[:, ::-1], -Jy)
plt.plot(x[ii[0]], y[ii[1]], "or", label="V = 1")
plt.title("Vector plot of current flow", fontsize=16)
plt.xlabel("X", fontsize=14)
plt.ylabel("Y", fontsize=14)
grid()
plt.legend()
# plt.savefig("7.png")

# Surface plot of potential
fig8 = figure(8, figsize=(7, 7))
ax = p3.Axes3D(fig8, auto_add_to_figure=False)
fig8.add_axes(ax)
plt.title("The 3-D surface plot of the potential", fontsize=16)
surface = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=cm.jet)
# plt.savefig("8.png")

plt.show()
