"""

Author: Jnaneswara Rao Rompilli
Date: 19-04-2022
Description: Spectra of non-periodic signals

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import *

# Estimage omega
def est_omega(w, Y):
    ii = np.where(w > 0)
    omega = sum(abs(Y[ii]) ** 2 * w[ii]) / sum(abs(Y[ii]) ** 2)  # weighted average
    return omega


# Estimate delta
def est_delta(w, Y):
    sup = 1e-4
    window = 1
    ii_1 = np.where(np.logical_and(np.abs(Y) > sup, w > 0))[0]
    np.sort(ii_1)
    points = ii_1[1 : window + 1]
    return np.sum(np.angle(Y[points])) / len(points)


# Function to plot spectrum after calculating
def plot_spectrum(index, w, Y, xlim, title, xlabel, ylabel1, ylabel2):
    plt.figure(index)
    plt.subplot(2, 1, 1)
    plt.plot(w, abs(Y), lw=2)
    plt.xlim([-xlim, xlim])
    plt.ylabel(ylabel1, fontsize=10)
    plt.title(title)
    plt.grid(True)
    plt.subplot(2, 1, 2)
    plt.plot(w, np.angle(Y), "ro", lw=2)
    plt.xlim([-xlim, xlim])
    plt.xlabel(xlabel, fontsize=10)
    plt.ylabel(ylabel2, fontsize=10)
    plt.grid(True)
    pass


# Intializing variables
PI = np.pi

# Q1: Spectrum of sin(sqrt(2)t), in the basic way
t = np.linspace(-PI, PI, 65)
t = t[:-1]
dt = t[1] - t[0]
fmax = 1 / dt
y = np.sin(np.sqrt(2) * t)
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y)) / 64.0
w = np.linspace(-PI * fmax, PI * fmax, 65)
w = w[:-1]


plot_spectrum(
    1,
    w,
    Y,
    10,
    r"Spectrum of $\sin(\sqrt{2}t)$",
    r"$\omega$",
    r"$|Y|$",
    r"Phase of $Y$",
)

# Spectrum of sin(sqrt(2)t), after windowing (better way)
t = np.linspace(-4 * PI, 4 * PI, 257)
t = t[:-1]
dt = t[1] - t[0]
n = np.arange(256)
wnd = fftshift(0.54 + 0.46 * np.cos(2 * PI * n / 256))
y = np.sin(np.sqrt(2) * t)
y = y * wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y)) / 256.0
w = np.linspace(-PI * fmax, PI * fmax, 257)
w = w[:-1]

plot_spectrum(
    2,
    w,
    Y,
    5,
    r"Spectrum of $\sin(\sqrt{2}t)$",
    r"$\omega$",
    r"$|Y|$",
    r"Phase of $Y$",
)

# Q2: f(t) = cos^3(w*t) with and without hamming window
t = np.linspace(-4 * PI, 4 * PI, 257)
t = t[:-1]
dt = t[1] - t[0]
fmax = 1 / dt
n = np.arange(256)
wnd = fftshift(0.54 + 0.46 * np.cos(2 * PI * n / 256))

# Without hamming window
y = np.cos(0.86 * t) * np.cos(0.86 * t) * np.cos(0.86 * t)
y[0] = 0
y1 = fftshift(y)
Y1 = fftshift(fft(y)) / 256.0

# With hamming window
y = np.cos(0.86 * t) * np.cos(0.86 * t) * np.cos(0.86 * t)
y = y * wnd
y[0] = 0
y2 = fftshift(y)
Y2 = fftshift(fft(y2)) / 256.0

w = np.linspace(-PI * fmax, PI * fmax, 257)
w = w[:-1]

plot_spectrum(
    3,
    w,
    Y1,
    10,
    r"Spectrum of $\cos^3(t)$ (without hamming window)",
    r"$\omega$",
    r"$|Y|$",
    r"Phase of $Y$",
)

plot_spectrum(
    4,
    w,
    Y2,
    7.5,
    r"Spectrum of $\cos^3(t)$ (with hamming window)",
    r"$\omega$",
    r"$|Y|$",
    r"Phase of $Y$",
)

# Q3: f(t) = cos(ωt + 𝛿)
# Here,
w0 = 1.2
d = 0.5

t = np.linspace(-PI, PI, 129)
t = t[:-1]
dt = t[1] - t[0]
fmax = 1 / dt
n = np.arange(128)
wnd = fftshift(0.54 + 0.46 * np.cos(2 * PI * n / 128))
y = np.cos(w0 * t + d)
y = y * wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y)) / 128.0

w = np.linspace(-PI * fmax, PI * fmax, 129)
w = w[:-1]

plot_spectrum(
    5, w, Y, 4, r"Spectrum of $\cos(ωt + 𝛿)$", r"$\omega$", r"$|Y|$", r"Phase of $Y$",
)

w_cal = est_omega(w, Y)
delta = est_delta(w, Y)
print("Calculated ω without noise: ", w_cal)
print("Calculated 𝛿 without noise: ", delta)

# Q4: After adding noise to the above question (Q3)
y = np.cos(w0 * t + d) + 0.1 * np.random.randn(128)
y = y * wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y)) / 128.0

plot_spectrum(
    6,
    w,
    Y,
    4,
    r"Spectrum of $\cos(ωt + 𝛿)$ (with noise)",
    r"$\omega$",
    r"$|Y|$",
    r"Phase of $Y$",
)

w_cal = est_omega(w, Y)
delta = est_delta(w, Y)
print("Calculated ω with noise: ", w_cal)
print("Calculated 𝛿 with noise: ", delta)

# Q5: Spectrum of chirped signal
t = np.linspace(-PI, PI, 1025)
t = t[:-1]
dt = t[1] - t[0]
fmax = 1 / dt
n = np.arange(1024)
wnd = fftshift(0.54 + 0.46 * np.cos(2 * PI * n / 1024))
y = np.cos(16.0 * (1.5 + t / (2.0 * PI)) * t)
y = y * wnd
y[0] = 0
y = fftshift(y)
Y = fftshift(fft(y)) / 1024.0

w = np.linspace(-PI * fmax, PI * fmax, 1025)
w = w[:-1]

plot_spectrum(
    7, w, Y, 100, r"Spectrum of chirped signal", r"$\omega$", r"$|Y|$", r"Phase of $Y$",
)

# Q6: Visualing spectrum evolution with time
t = np.linspace(-PI, PI, 1025)
t = t[:-1]
dt = t[1] - t[0]
fmax = 1 / dt

# Splitting time array
t_array = np.split(t, 16)
Y_mag = np.zeros((16, 64))
Y_phase = np.zeros((16, 64))

# Calculating DFT for each piece
for i in range(len(t_array)):
    n = np.arange(64)
    wnd = fftshift(0.54 + 0.46 * np.cos(2 * PI * n / 64))
    y = np.cos(16.0 * (1.5 + t_array[i] / (2.0 * PI)) * t_array[i])
    y = y * wnd
    y[0] = 0
    y = fftshift(y)
    Y = fftshift(fft(y)) / 64.0

    Y_mag[i] = abs(Y)
    Y_phase[i] = np.angle(Y)

t = t[::64]
w = np.linspace(-fmax * PI, fmax * PI, 65)
w = w[:-1]
t, w = np.meshgrid(t, w)

# 3D Plots - Magnitude, Phase
fig = plt.figure(8)
plt.title("Time-frequency plot")
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(w, t, Y_mag.T, cmap="viridis", linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.ylabel("Frequency")
plt.xlabel("time")

fig = plt.figure(9)
plt.title("Time-frequency plot")
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(w, t, Y_phase.T, cmap="viridis", linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.ylabel("Frequency")
plt.xlabel("time")


plt.show()
