
# The Flow Past an Oscillating Circular Cylinder in a Channel

$$
|\mathbf{U}\_\infty| = 2 \pi f A_{\text{ref}}
                     = 2 \times \pi \times 0.25 \times 0.25
               \approx 0.3927\ \frac{\mathrm{m}}{\mathrm{s}}
$$
$$
p\_{\text{dyn}} = \frac{1}{2} \rho_{\infty} \mathbf{U}\_\infty^2
                = 0.5 \times 1.0 \times (0.3927)^2
          \approx 0.0771\ \mathrm{Pa}
$$
$$
A_{\text{ref}} = \Delta z\ l_{\text{ref}} = 0.1 \times 0.1 = 0.01 \mathrm{m}^2
$$
$$
\texttt{forceScaling} = \frac{1}{A_{\text{ref}} p_{\text{dyn}}}
                      = \frac{1}{0.01 \times 0.0771}
                \approx 1297\ \mathrm{N}^{-1}
$$

This means that the total fluid forces (in the $x$- and $y$-directions) acting
on the cylinder are scaled by a factor with units of $\mathrm{N}^{-1}$, in this
case, $1927\ \mathrm{N}^{-1}$, to ultimately yield the force coefficients.

$$
C_{\mathrm{d}} = \texttt{forceScaling}\ F_x
$$
$$
C_{\mathrm{l}} = \texttt{forceScaling}\ F_y
$$
