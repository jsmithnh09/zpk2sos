# zpk2sos
A Zero-Pole-Gain to Second-Order-Section conversion library.

## Background
With most filter design applications, the poles and zeros can be generated with different design parameters, resulting in a zero/pole/gain filter structure. In real world applications, these should be implemented as Biquads that have dedicated hardware routines.

## Algorithm
The algorithm tackles a couple primary steps:
1. Sort the poles such that each singularity is adjacent to its conjugate
3. Match poles closest to unit circle, `|x|+eps >= 1`, with zeros closest to those poles.
4. Re-iterate #3 until every singularity has a paired pole or zero.
5. For real poles, group them together (per section) that are closest in absolute value. The same rule holds for real zeros, (reduces the number of large/small value multiplications that result in roundoff error.)

## References
\[1\] Simon Kornblith, Galen Lynch, Martin Holters, João Felipe Santos, Spencer Russell, Jay Kickliter, Jeff Bezanson, Gudmundur Adalsteinsson, Alex Arslan, Ryuichi Yamamoto, jordancluts, Matti Pastell, Tony Kelman, Ben Arthur, Tom Krauss, HDictus, Hamza El-Saawy, Jared Kofron, Eric Hanson, … Jordan Smith. (2023). JuliaDSP/DSP.jl: v0.7.9 (v0.7.9). Zenodo. https://doi.org/10.5281/zenodo.8344531
\[2\] Virtanen, P., Gommers, R., Oliphant, T.E. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020). https://doi.org/10.1038/s41592-019-0686-2
