"""
Author: Jnaneswara Rao Rompilli
Date: 15-04-2022
Description: The Digital Fourier Transform

"""
# Importing necesaary modules
from scipy.fft import *
import numpy as np
import matplotlib.pyplot as plt

PI = np.pi

# FFT of random values
x = np.random.rand(100)
X = fft(x)
y = ifft(X)
np.c_[x, y]
print(abs(x - y).max())

# f(t) = sin(5*t)
x = np.linspace(0, 2 * PI, 129)
x = x[:-1]
y = np.sin(5 * x)
Y = fftshift(fft(y)) / 128.0
w = np.linspace(-64, 63, 128)

# Magnitude and phase plot of above function
plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(w, abs(Y), lw=2)
plt.xlim([-10, 10])
plt.ylabel(r"$|Y|$", fontsize=10)
plt.title(r"Spectrum of $\sin(5t)$", fontsize=12)
plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(w, np.angle(Y), "ro", lw=2)
ii = np.where(abs(Y) > 1e-3)
plt.plot(w[ii], np.angle(Y[ii]), "go", lw=2)
plt.xlim(-10, 10)
plt.ylabel(r"Phase of $Y$", fontsize=10)
plt.xlabel(r"$k$", fontsize=10)
plt.grid(True)


# f(t) = (1 + 0.1cos(t))*cos(10t)
t = np.linspace(-4 * PI, 4 * PI, 513)
t = t[:-1]
y = (1 + 0.1 * np.cos(t)) * np.cos(10 * t)
Y = fftshift(fft(y)) / 512.0
w = np.linspace(-64, 63, 513)
w = w[:-1]

# Magnitude and phase plot of above function
plt.figure(2)
plt.subplot(2, 1, 1)
plt.plot(w, abs(Y), lw=2)
plt.xlim([-15, 15])
plt.ylabel(r"|Y|", fontsize=10)
plt.title(
    r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$",
    fontsize=12,
)
plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(w, np.angle(Y), "ro", lw=2)
plt.xlim([-15, 15])
plt.ylabel(r"Phase of $Y$", fontsize=10)
plt.xlabel(r"$\omega$", fontsize=10)
plt.grid(True)


# f(t) = sin^3(t) = (3*sin(t) - sin(3*t)) / 4
t = np.linspace(-4 * PI, 4 * PI, 513)
t = t[:-1]
y = np.sin(t) * np.sin(t) * np.sin(t)
Y = fftshift(fft(y)) / 512.0
w = np.linspace(-64, 64, 513)
w = w[:-1]

# Magnitude plot of above function
plt.figure(3)
plt.plot(w, abs(Y), lw=2)
plt.xlim([-10, 10])
plt.xlabel(r"$\omega$", fontsize=10)
plt.ylabel(r"|Y|", fontsize=10)
plt.title("Spectrum of sin^3(t)", fontsize=12)
plt.grid(True)


# f(t) = cos^3(t) = (3*cos(t) + cos(3*t)) / 4
t = np.linspace(-4 * PI, 4 * PI, 513)
t = t[:-1]
y = np.cos(t) * np.cos(t) * np.cos(t)
Y = fftshift(fft(y)) / 512.0
w = np.linspace(-64, 64, 513)
w = w[:-1]

# Magnitude plot of above function
plt.figure(4)
plt.plot(w, abs(Y), lw=2)
plt.xlim([-10, 10])
plt.ylabel(r"|Y|", fontsize=10)
plt.xlabel(r"$\omega$", fontsize=10)
plt.title("Spectrum of cos^3(t)", fontsize=12)
plt.grid(True)


# f(t) = cos(20t + 5cos(t))
t = np.linspace(-4 * PI, 4 * PI, 513)
t = t[:-1]
y = np.cos(20 * t + 5 * np.cos(t))
Y = fftshift(fft(y)) / 512.0
w = np.linspace(-64, 64, 513)
w = w[:-1]

# Magnitude and phase plot of above function
plt.figure(5)
plt.subplot(2, 1, 1)
plt.plot(w, abs(Y), lw=2)
plt.xlim([-40, 40])
plt.ylabel(r"$|Y|$", fontsize=10)
plt.title(r"Spectrum of $\cos(20t + 5cos(t))$", fontsize=12)
plt.grid(True)
plt.subplot(2, 1, 2)
ii = np.where(abs(Y) > 1e-3)
plt.plot(w[ii], np.angle(Y[ii]), "go", lw=2)
plt.xlim(-40, 40)
plt.ylabel(r"Phase of $Y$", fontsize=10)
plt.xlabel(r"$\omega$", fontsize=10)
plt.grid(True)


# f(t) = exp(-t*t/2)
N = 128
T = 2 * PI
tolerance = 1e-15
error = 1
w = []
Y = []
# Only when the error is less than a tolerance value will the loop be terminated
while True:
    t = np.linspace(-T / 2, T / 2, N + 1)
    t = t[:-1]
    w = N / T * np.linspace(-PI, PI, N + 1)
    w = w[:-1]

    y = np.exp(-0.5 * t * t)
    Y = fftshift(fft(y)) * T / (2 * PI * N)
    Y_act = (1 / np.sqrt(2 * PI)) * np.exp(-0.5 * w * w)

    error = np.mean(abs(abs(Y) - Y_act))

    if error < tolerance:
        break

    T, N = 2 * T, 2 * N

print(f"For accurate frequency domain: N = {N}, T = {T/PI}*π")
# Magnitude plot of above function
plt.figure(6)
plt.plot(w, abs(Y), lw=2)
plt.xlim([-10, 10])
plt.ylabel(r"$|Y|$", fontsize=10)
plt.xlabel(r"$\omega$", fontsize=10)
plt.title("Spectrum of Gaussian function", fontsize=12)
plt.grid(True)


plt.show()
