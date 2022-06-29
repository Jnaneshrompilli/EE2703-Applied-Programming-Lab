"""-------------------------------------------------------------------

Author: Jnaneswara Rao Rompilli (EE20B052)
Date: 22-2-22, 22:22
Description: Calculating fourier coefficients of exp(x) & cos(cos(x))

---------------------------------------------------------------------"""

from pylab import *
from scipy import integrate

# Declaring constants
STEP_SIZE = 0.1
PI = np.pi

# Generating samples in a range for plotting
def generateSample(x1, x2):
    sample = []

    while x1 <= x2:
        sample.append(x1)
        x1 += STEP_SIZE

    return np.array(sample)


# Function for exp(x)
def expo(t):
    return exp(t)


# Function for cos(cos(x))
def coscos(t):
    return cos(cos(t))


# Functions to pass in integrate.quad function for calculating Fourier coefficient values
def a_coscos(t, k):
    return cos(cos(t)) * cos(k * t)


def b_coscos(t, k):
    return cos(cos(t)) * sin(k * t)


def a_expo(t, k):
    return exp(t) * cos(k * t)


def b_expo(t, k):
    return exp(t) * sin(k * t)


# To Plot both functions in given range
x_plot = generateSample(-2 * PI, 4 * PI)

# For fitting the coefficients in this range
x = linspace(0, 2 * PI, 401)
x = x[:-1]

# Computing function values
y_exp = expo(x)
y_coscos = coscos(x)

# To store fourier coefficients calculated from integration
coeff_exp = np.zeros((51,1))
coeff_coscos = np.zeros((51,1))

an_exp, bn_exp, an_coscos, bn_coscos = 0, 0, 0, 0
index = 0

# Computing a0 seperately
an_exp = integrate.quad(a_expo, 0.0, 2.0 * PI, args=(0))[0] / (2 * PI)
an_coscos = integrate.quad(a_coscos, 0.0, 2.0 * PI, args=(0))[0] / (2 * PI)

coeff_exp[index] = an_exp
coeff_coscos[index] = an_coscos
index += 1

n = np.arange(1, 52)

# Calculating each coefficient(a1,b1,a2,b2,...) value using integrate.quad() function
for k in range(1, 26):
    coeff_exp[index] = integrate.quad(a_expo, 0.0, 2.0 * PI, args=(k))[0] / PI
    coeff_exp[index + 1] = integrate.quad(b_expo, 0.0, 2.0 * PI, args=(k))[0] / PI

    coeff_coscos[index] = integrate.quad(a_coscos, 0.0, 2.0 * PI, args=(k))[0] / PI
    coeff_coscos[index + 1] = integrate.quad(b_coscos, 0.0, 2.0 * PI, args=(k))[0] / PI

    index += 2


# Fitting values and calculating fourier coefficients
# x = x[:-1]
count = len(x)
M = np.zeros((count, 51))
M[:, 0] = 1
for k in range(1, 26):
    index = 2 * k
    M[:, index - 1] = cos(k * x)
    M[:, index] = sin(k * x)

coefflst_coscos = lstsq(M, y_coscos, rcond=None)[0].reshape((-1, 1))
coefflst_exp = lstsq(M, y_exp, rcond=None)[0].reshape((-1, 1))

# Finding deviation from actual values
dev_exp = abs(coefflst_exp - coeff_exp).max()
dev_coscos = abs(coefflst_coscos - coeff_coscos).max()

print("Maximum deviation for exp(x): ", dev_exp)
print("Maximum deviation for cos(cos(x)): ", dev_coscos)

# Computing y values from calculated fourier coefficients
ylst_exp = np.dot(M, coefflst_exp)
ylst_coscos = np.dot(M, coefflst_coscos)

# Plotting actual value and values calculated from fourier coefficients for exp(x)
plt.figure(1, figsize=(8, 6))
plt.semilogy(x_plot, expo(x_plot), color="darkorange", label="Actual", linewidth=3)
plt.semilogy(x[::5], ylst_exp[::5], "go", markersize=4, label="Fourier prediction")
plt.title("Plotting exponential exp(x)", size=20)
plt.xlabel(f"x", size=20)
plt.ylabel(f"f(x)", size=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend(loc="best")

# Plotting actual value and values calculated from fourier coefficients for cos(cos(x))
plt.figure(2, figsize=(8, 6))
plt.plot(x_plot, coscos(x_plot), color="darkorange", label="Actual", linewidth=3)
plt.plot(x[::5], ylst_coscos[::5], "go", markersize=4, label="Fourier prediction")
plt.title("Plotting cosine", size=20)
plt.xlabel(f"x", size=20)
plt.ylabel(f"f(x)", size=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend(loc="best")

# Plotting fourier coefficients for exp(x) [actual and predicted] in semilogy scale
plt.figure(3, figsize=(8, 6))
plt.semilogy(n, abs(coeff_exp), "or", label="Integration")
plt.semilogy(n, abs(coefflst_exp), "og", label="Prediction")
plt.title("Fourier coefficients of exp(x)", size=20)
plt.xlabel(f"n", size=20)
plt.ylabel(f"Coefficients(log)", size=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend(loc="best")

# Plotting fourier coefficients for exp(x) [actual and predicted] in loglog scale
plt.figure(4, figsize=(8, 6))
plt.loglog(n, abs(coeff_exp), "or", label="Integration")
plt.loglog(n, abs(coefflst_exp), "og", label="Prediction")
plt.title("Fourier coefficients of exp(x)", size=20)
plt.xlabel(f"n (log)", size=20)
plt.ylabel(f"Coefficients(log)", size=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend(loc="best")

# Plotting fourier coefficients for cos(cos(x)) [actual and predicted] in semilogy scale
plt.figure(5, figsize=(8, 6))
plt.semilogy(n, abs(coeff_coscos), "or", label="Integration")
plt.semilogy(n, abs(coefflst_coscos), "og", label="Prediction")
plt.title("Fourier coefficients of cos(cos(x))", size=20)
plt.xlabel(f"n", size=20)
plt.ylabel(f"Coefficients(log)", size=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend(loc="best")

# Plotting fourier coefficients for cos(cos(x)) [actual and predicted] in loglog scale
plt.figure(6, figsize=(8, 6))
plt.loglog(n, abs(coeff_coscos), "or", label="Integration")
plt.loglog(n, abs(coefflst_coscos), "og", label="Prediction")
plt.title("Fourier coefficients of cos(cos(x))", size=20)
plt.xlabel(f"n (log)", size=20)
plt.ylabel(f"Coefficients(log)", size=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True)
plt.legend(loc="best")

plt.show()
