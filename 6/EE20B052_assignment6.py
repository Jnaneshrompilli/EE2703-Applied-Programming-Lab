""" ******************************************************
Author: Jnaneswara Rao Rompilli (EE20B052)
Date: 26-03-2022
Description: Analyze Linear Time-Invariant systems using Laplace Transform

********************************************************* """

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sp

# Plot impulse response for given transfer function with different decay values
def plot_impulse(decay, level, cnt):
    H = sp.lti(
        [1, decay], np.polymul([1, 0, 2.25], [1, 2 * decay, 2.25 + decay * decay])
    )  # Transfer function in s-domain

    w, S, phi = H.bode()  # Bode
    plt.figure(cnt + 1)
    plt.subplot(2, 1, 1)  #  Magnitude plot
    plt.title(f"Bode plot of transfer function ({level} decay)", fontsize=12)
    plt.semilogx(w, S, label="|H(s)|")
    plt.xlabel("ω (log scale)", fontsize=10)
    plt.ylabel("|H(s)| (decibel)")
    plt.legend()
    plt.grid()
    plt.subplot(2, 1, 2)  # Pase plot
    plt.semilogx(w, phi, label="<H(jw)")
    plt.xlabel("ω (log scale)", fontsize=10)
    plt.ylabel("<H(s)")
    plt.legend()
    plt.grid()
    plt.savefig(f"{cnt+1}.png")

    # Impulse response
    t, x = sp.impulse(H, None, np.linspace(0, 50, 500))
    plt.figure(cnt + 2)
    plt.plot(t, x)
    plt.title(f"Time response of spring (decay = {decay})")
    plt.xlabel("t", fontsize=10)
    plt.ylabel("x(t)", fontsize=10)
    plt.grid()
    plt.savefig(f"{cnt+2}.png")


# 1, 2: Plotting impulse response for two different decay values
plot_impulse(0.5, "High", 0)
plot_impulse(0.05, "Low", 2)

# 3: Plotting responses for various frequency values
H = sp.lti([1], [1, 0, 2.25])  # Transfer function in s-domain
omega = 1.4
while omega <= 1.6:
    t = np.linspace(0, 50, 500)
    f_t = np.cos(omega * t) * np.exp(-0.05 * t)  # Different omega values
    t, x, svec = sp.lsim(H, f_t, t)

    # Output response
    plt.figure(5)
    plt.plot(t, x, label=f"w = {omega}")
    plt.title("x(t) for in frequency range", fontsize=12)
    plt.xlabel("t", fontsize=10)
    plt.ylabel("x(t)", fontsize=10)
    plt.legend()
    plt.grid()
    plt.savefig("5.png")

    omega += 0.05

# 4: Coupled spring problem
X = sp.lti([1, 0, 2], [1, 0, 3, 0])  # Transfer fucntions in s-domain
Y = sp.lti([2], [1, 0, 3, 0])

# Calculating x(t) and y(t) solutions for spring problem
plt.figure(6)
t, x = sp.impulse(X, None, np.linspace(0, 20, 500))
plt.plot(t, x, label="x(t)")
t, y = sp.impulse(Y, None, np.linspace(0, 20, 500))
plt.plot(t, y, label="y(t)")
plt.title(f"Coupled spring solution", fontsize=12)
plt.xlabel("t", fontsize=10)
plt.ylabel("function", fontsize=10)
plt.legend(loc="upper right")
plt.grid()
plt.savefig("6.png")

# 5: Two port network
# Values
L = 1e-6
C = 1e-6
R = 100

H = sp.lti([1], [L * C, R * C, 1])  # Transfer function in s-domain
w, S, phi = H.bode()
plt.figure(7)
plt.subplot(2, 1, 1)  #  Magnitude
plt.title(f"Steady State Transfer function of two-port network)", fontsize=12)
plt.semilogx(w, S, label="|H(s)|")
plt.xlabel("ω (log scale)", fontsize=10)
plt.ylabel("|H(s)| (decibel)")
plt.legend()
plt.grid()
plt.subplot(2, 1, 2)  # Angle
plt.semilogx(w, phi, label="<H(jw)")
plt.xlabel("ω (log scale)", fontsize=10)
plt.ylabel("<H(s)")
plt.legend()
plt.grid()
plt.savefig("7.png")

# 6: Output of Vi(t)
t = np.arange(0, 20e-3, 1e-7)
v_t = np.cos(1e3 * t) - np.cos(1e6 * t)
t, x, svec = sp.lsim(H, v_t, t)

# Long time scale plot of Output
plt.figure(8)
plt.plot(t, x)
plt.title("Output voltage of RLC two-port network (long time)", fontsize=12)
plt.xlabel("t", fontsize=10)
plt.ylabel("Vo(t)", fontsize=10)
plt.grid()
plt.savefig("8.png")

# Short time scale plot of Output
plt.figure(9)
plt.plot(t[0:500], x[:500])
plt.title("Output voltage of RLC two-port network (short time)", fontsize=12)
plt.xlabel("t", fontsize=10)
plt.ylabel("Vo(t)", fontsize=10)
plt.grid()
plt.savefig("9.png")


# plt.show()
