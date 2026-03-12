#!/usr/bin/python3
import numpy as np
from scipy import signal

# example: 8th‑order elliptic low‑pass, 1 dB ripple, 40 dB stop, Wn=0.2
N, rp, rs, Wn = 8, 1, 40, 0.2
z, p, k = signal.ellip(N, rp, rs, Wn, output='zpk')
sos = signal.zpk2sos(z, p, k, pairing='nearest')

np.set_printoptions(precision=10, suppress=True)
print(sos)
