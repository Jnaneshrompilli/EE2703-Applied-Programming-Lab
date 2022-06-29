"""-------------------------------------------------------

Author: Jnaneswara Rao Rompilli (EE20B052)
Date: 08-04-2022
Description: Analysis of circuits using Laplace Transforms

----------------------------------------------------------"""

from __future__ import division
import scipy.signal as sp
import numpy as np
from sympy import *
import pylab

PI = np.pi
s = symbols("s")

# Calculate frequency response (Vo) of Low pass filter
def lowpass(R1, R2, C1, C2, G, Vi):
    A = Matrix(
        [
            [0, 0, 1, -1 / G],
            [-1 / (1 + s * R2 * C2), 1, 0, 0],
            [0, -G, G, 1],
            [-1 / R1 - 1 / R2 - s * C1, 1 / R2, 0, s * C1],
        ]
    )
    b = Matrix([0, 0, 0, -Vi / R1])

    V = A.inv() * b  # [V1 Vp Vm Vo]

    return (A, b, V)


# Calculate frequency response (Vo) of High pass filter
def highpass(R1, R3, C1, C2, G, Vi):
    A = Matrix(
        [
            [0, -1, 0, 1 / G],
            [s * C2 * R3 / (s * C2 * R3 + 1), 0, -1, 0],
            [0, G, -G, 1],
            [-s * C2 - 1 / R1 - s * C1, 0, s * C2, 1 / R1],
        ]
    )

    b = Matrix([0, 0, 0, -Vi * s * C1])

    V = A.inv() * b

    return (A, b, V)


# Convert symbolic expression individual coefficients
def sympytolti(Hs):
    # Hs = expand(simplify(Hs))
    n, d = fraction(Hs)
    num = Poly(n, s)
    den = Poly(d, s)

    num_c = num.all_coeffs()
    den_c = den.all_coeffs()

    num_c, den_c = [float(f) for f in num_c], [float(f) for f in den_c]

    H = sp.lti(num_c, den_c)

    return H


# Calculate response of High Pass Filter to sinusoid
def hp_sinusoid(b):
    A3, b3, V3 = highpass(10000, 10000, 1e-9, 1e-9, 1.586, 1)
    Vo = V3[3]
    H3 = sympytolti(Vo)

    a = -500
    b = b * PI

    t = np.linspace(0, 1e-2, 100000)
    Vin = np.exp(a * t) * np.cos(b * t)
    t, Vt, svec = sp.lsim(H3, Vin, t)

    return (t, Vt, svec)


# Calculate response of Low Pass Filter to sinusoid
def lp_sinusoid(b):
    A4, b4, V4 = lowpass(10000, 10000, 1e-9, 1e-9, 1.586, 1)
    Vo = V4[3]
    H4 = sympytolti(Vo)

    a = -500
    b = b * PI

    t = np.linspace(0, 1e-2, 100000)
    Vin = np.exp(a * t) * np.cos(b * t)
    t, Vt, svec = sp.lsim(H4, Vin, t)

    return (t, Vt, svec)


# Low pass filter
A, b, V = lowpass(10000, 10000, 1e-9, 1e-9, 1.586, 1)
Vo = V[3]
H = sympytolti(Vo)

ww = pylab.logspace(0, 8, 801)
ss = 1j * ww
hf = lambdify(s, Vo, "numpy")
v = hf(ss)

pylab.figure(1, figsize=(7, 6))
pylab.loglog(ww, abs(v), lw=2)
pylab.title("Low Pass filter", fontsize=12)
pylab.xlabel("ω")
pylab.ylabel("Magnitude")
pylab.grid(True)

# Step response of a Low pass filter
A1, b1, H1 = lowpass(10000, 10000, 1e-9, 1e-9, 1.586, 1 / s)
Vo = H1[3]
Vh = sympytolti(Vo)

t, Vt = sp.impulse(Vh, None, np.linspace(0, 5e-3, 10000))
pylab.figure(2, figsize=(7, 6))
pylab.plot(t, Vt, label="V(t)")
pylab.title("Step response of Low Pass filter", fontsize=12)
pylab.xlabel("t", fontsize=10)
pylab.ylabel("Vo(t)", fontsize=10)
pylab.grid()
pylab.legend()



# Output voltage for given Vi(t)
Vit = np.sin(2000 * PI * t) + np.cos(2e6 * PI * t)
t, Vo, svec = sp.lsim(H, Vit, t)

pylab.figure(3, figsize=(7, 6))
pylab.plot(t, Vo)
pylab.title("Output voltage Vo(t) (Low pass)", fontsize=12)
pylab.xlabel("t")
pylab.ylabel("Vo(t)")
pylab.grid(True)



# High Pass filter
A2, b2, H2 = highpass(10e3, 10e3, 1e-9, 1e-9, 1.586, 1)
Vo = H2[3]

ww = pylab.logspace(0, 8, 801)
ss = 1j * ww
hf = lambdify(s, Vo, "numpy")
v = hf(ss)

pylab.figure(4, figsize=(7, 6))
pylab.loglog(ww, abs(v), lw=2)
pylab.title("High Pass filter", fontsize=12)
pylab.xlabel("ω")
pylab.ylabel("Magnitude")
pylab.grid(True)


# Output response of a High pass filter for a damped sinusoid - Low frequency
t, Vt, svec = hp_sinusoid(2000)

pylab.figure(5, figsize=(7, 6))
pylab.plot(t, Vt, label="V(t)")
pylab.title(
    "Response of High Pass filter to damped sinusoid (Low frequency)", fontsize=12
)
pylab.xlabel("t", fontsize=10)
pylab.ylabel("Vo(t)", fontsize=10)
pylab.grid()
pylab.legend()

# Output response of a Low pass filter for a damped sinusoid - Low frequency
t, Vt, svec = lp_sinusoid(2000)
pylab.figure(6, figsize=(7, 6))
pylab.plot(t, Vt, label="V(t)")
pylab.title(
    "Response of Low Pass filter to damped sinusoid (Low frequency)", fontsize=12
)
pylab.xlabel("t", fontsize=10)
pylab.ylabel("Vo(t)", fontsize=10)
pylab.grid()
pylab.legend()

# Output response of a High pass filter for a damped sinusoid - High frequency
t, Vt, svec = hp_sinusoid(2000000)
pylab.figure(7, figsize=(7, 6))
pylab.plot(t, Vt, label="V(t)")
pylab.title(
    "Response of High Pass filter to damped sinusoid (High frequency)", fontsize=12
)
pylab.xlabel("t", fontsize=10)
pylab.ylabel("Vo(t)", fontsize=10)
pylab.grid()
pylab.legend()

# Output response of a Low pass filter for a damped sinusoid - High frequency
t, Vt, svec = lp_sinusoid(2000000)
pylab.figure(8, figsize=(7, 6))
pylab.plot(t, Vt, label="V(t)")
pylab.title(
    "Response of Low Pass filter to damped sinusoid (High frequency)", fontsize=12
)
pylab.xlabel("t", fontsize=10)
pylab.ylabel("Vo(t)", fontsize=10)
pylab.grid()
pylab.legend()



# Step response of a High pass filter
A3, b3, H3 = highpass(10000, 10000, 1e-9, 1e-9, 1.586, 1 / s)
Vo = H3[3]
Vh = sympytolti(Vo)

t, Vt = sp.impulse(Vh, None, np.linspace(0, 5e-3, 10000))
pylab.figure(9, figsize=(7, 6))
pylab.plot(t, Vt, label="V(t)")
pylab.title("Step response of High pass filter", fontsize=12)
pylab.xlabel("t", fontsize=10)
pylab.ylabel("Vo(t)", fontsize=10)
pylab.grid()
pylab.legend()

pylab.show()

