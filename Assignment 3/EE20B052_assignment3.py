from pylab import *
import scipy.special as sp


# Intializing coefficients
A = 1.05
B = -0.105
k = 9

# Actual function
def g(t, A, B):
    value = A * sp.jn(2, t) + B * t
    return value


# Calculate RMS error of A and B values w.r.t to Y1
def rms_error(data, t, a, b):
    n = len(t)
    error_value = 0
    for k in range(0, n):
        error_value += ((data[:, 1][k] - g(t[k], a, b)) ** 2) / n
    return error_value


data = []

# Take input from .dat file
try:
    data = loadtxt("fitting.dat")
except IOError:
    print("File Path doesn't exist !")
    exit()


t = np.array(data[:, 0])  # X co-ordinates
y_actual = g(t, A, B)
stdev = logspace(-1, -3, k)  # Calculate standard deviation
n, m = data.shape


# Plot actual function and noise
figure(0, figsize=(8, 6))
for i in range(0, k):
    plot(t, data[:, i + 1], label=f"σ={round(stdev[i],3)}")
plot(t, y_actual, color="black", label="Actual function", linewidth=3)
title("Data Plot", size=20)
xlabel(f"t", size=20)
ylabel(f"f(t) + n", size=20)
xticks(fontsize=12)
yticks(fontsize=12)
grid(True)
legend(loc="best")
# plt.savefig("0.png")

# Plot first column with error bars
figure(1, figsize=(8, 6))
errorbar(t[::5], data[:, 1][::5], stdev[0], fmt="ro")
plot(t, y_actual, color="black", label="Actual function", linewidth=3)
title("Data with Error Bars", size=20)
xlabel(f"t", size=20)
ylabel(f"f(t) + n(0, {stdev[0]})", size=20)
xticks(fontsize=12)
yticks(fontsize=12)
grid(True)
legend(loc="best")
# plt.savefig("1.png")

# Calculate M matrix
M = np.array((n, 2))
M = c_[sp.jn(2, t), t]
p = [A, B]
y = dot(M, p)

if array_equal(y, y_actual):
    print("Both column vectors are equal")
else:
    print("Both column vectors are not equal")

# Calculate RMS error of A, B value in a range
A1 = linspace(0, 2, 21, endpoint=True)
B1 = linspace(-0.2, 0, 21, endpoint=True)
E = np.zeros((21, 21))

# Find rms error of Given A and B values
for i in range(0, 21):
    for j in range(0, 21):
        E[i][j] += rms_error(data, t, A1[i], B1[j])

X, Y = np.meshgrid(A1, B1)

# Plot contours of error
figure(2, figsize=(8, 6))
ctr = contour(X, Y, E, 20)
clabel(ctr, ctr.levels[0:5], inline=True, fontsize=10)
title("Contour Plot", size=20)
xlabel("A", size=20)
ylabel("B", size=20)
xticks(fontsize=12)
yticks(fontsize=12)
grid(True)
# plt.savefig("2.png")


# Estimate A, B for each coloumn
prediction = np.zeros((k, 2))
error = np.zeros((k, 2))
for i in range(1, 10):
    AB = linalg.lstsq(M, data[:, i], rcond=None)
    prediction[i - 1] = AB[0]
    error[i - 1] = abs(prediction[i - 1] - [A, B])

# Plot estimation errors of A and B w.r.t standard deviation in normal scale
figure(3, figsize=(8, 6))
plot(stdev, error[:, 0], marker="o", label="A", linestyle="dashed")
plot(stdev, error[:, 1], marker="o", label="B", linestyle="dashed")
title("Error in predicted coefficients", size=20)
xlabel(f"Standard Deviation of noise", size=20)
ylabel(f"RMS error", size=20)
xticks(fontsize=12)
yticks(fontsize=12)
grid(True)
legend(loc="best")
# plt.savefig("3.png")

# Plot estimation errors of A and B w.r.t standard deviation in log scale
figure(4, figsize=(8, 6))
errorbar(stdev, error[:, 0], stdev, fmt="ro")
loglog(stdev, error[:, 0], marker="o", label="A", linestyle="none")
loglog(stdev, error[:, 1], marker="o", label="B", linestyle="none")
errorbar(stdev, error[:, 1], stdev, fmt="go")
title("Error in predicted coefficients", size=20)
xlabel(f"Standard Deviation of noise(log)", size=20)
ylabel(f"RMS error(log)", size=20)
xticks(fontsize=12)
yticks(fontsize=12)
grid(True)
legend(loc="best")
# plt.savefig("4.png")

show()

