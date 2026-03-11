# Quintom Model and the Shooting Method in CLASS-QTM

This document provides a detailed explanation of the Quintom model and its implementation using a shooting method, as seen in `modeling/model_good2.ipynb`.

## 1. The Quintom Model

The **Quintom model** is a dark energy scenario involving two scalar fields to provide a dynamic equation of state $w$ that can cross the phantom divide ($w = -1$).

### Scalar Fields
1.  **Quintessence field ($\phi$):** A canonical scalar field with a positive kinetic term.
2.  **Phantom field ($\sigma$):** A non-canonical scalar field with a negative kinetic term.

### Potential and Coupling
- **Potential:** The quintessence field $\phi$ typically follows an exponential potential:
  $$V(\phi) = V_0 e^{-\lambda_\phi \phi}$$
- **Coupling:** In this specific implementation, the phantom field $\sigma$ is coupled to the matter density $\rho_m$ (Baryons + CDM):
  $$\rho_m = K_m e^{\delta \sigma - 3N}$$
  where $\delta$ is the coupling constant and $N = \ln a$ is the number of e-folds.

### Analytical Calculation of $\Omega_\phi$ and $\Omega_\sigma$ (RD Era)

In the radiation-dominated (RD) era, the scalar fields follow scaling solutions. We can calculate their energy densities relative to the radiation density $\rho_{rad}$.

#### 1. Quintessence Field ($\phi$)
The energy density $\Omega_\phi$ in the RD era is determined by the potential $V(\phi)$ and the kinetic term. 
In the scaling regime ($a \ll a_{eq}$):
$$V(\phi) \propto \rho_{rad}$$
The analytical expression used in the notebook is:
$$V_{RD} = \frac{1}{3} Q(\delta, \lambda_\phi) \rho_{rad}$$
The total energy density for the quintessence field is $\rho_\phi = 3 V$, leading to the analytical ratio:
$$\Omega_\phi^{RD} \approx \frac{\rho_\phi^{RD}}{\rho_{rad}} = Q(\delta, \lambda_\phi)$$

#### 2. Phantom Field ($\sigma$)
The phantom field in this model has a strictly negative kinetic term and is coupled to matter. In the early universe, it is typically assumed to start from rest ($\sigma' = 0$), so:
$$\Omega_\sigma^{RD} \approx 0$$

#### 3. Total Dark Energy
The total dark energy density parameter is the sum of both components:
$$\Omega_{DE} = \Omega_\phi + \Omega_\sigma$$

---

## 2. The Shooting Method

The shooting method is a numerical technique used to solve boundary value problems. In cosmology, we know the desired parameters **today** (e.g., $\Omega_m, \Omega_{DE}$), but we must integrate the equations from the **early universe** ($a \approx 10^{-14}$).

### Why it is needed?
The initial values of $\phi$ and $\sigma$ in the radiation era are not direct observables. We must "shoot" (vary) the initial parameters or global constants until the late-time evolution matches our target observations.

### Shooting Variables (Unknowns)
These are the input parameters varied by the solver:
1.  **`logV0`**: The logarithm of the potential amplitude $V_0$. It primarily determines the energy density of the quintessence field today.
2.  **`Km`**: The matter scaling constant. Since $\rho_m \approx K_m$ today (if $\sigma$ is small), this is tuned to match the observed matter density $\Omega_{m,0}$.
3.  **`h`**: The dimensionless Hubble parameter. It is tuned to ensure the angular size of the sound horizon matches CMB data.

### Targets (Constraints)
These are the physical values the system must converge to:
1.  **`Flatness`**: Ensures that the total energy density $\sum \rho_i$ equals the critical density $H^2$ at $z=0$ (i.e., $\Omega_{tot} = 1$).
2.  **`Omega_de` (or `Omega_m`)**: The target dark energy parameter today (typically $\sim 0.69$).
3.  **`100*theta_s`**: The angular size of the sound horizon $\theta_s = r_s / D_A$. This is a high-precision target from Planck data.

---

## 3. How the Shooting Method Works

The process is managed by the `ShootingManager` class using the `scipy.optimize.least_squares` algorithm:

1.  **Initialization**: The solver starts with an initial guess for `logV0`, `Km`, and `h`.
2.  **Update**: It defines new cosmological parameters based on the current guess.
3.  **Simulation**:
    -   It calculates the **decoupling redshift** ($z_{dec}$) using a Peebles recombination model.
    -   It integrates the system of ODEs ($\phi, \phi', \sigma, \sigma'$) from $N_{ini}$ to $N=0$.
4.  **Derived Parameters**:
    -   It calculates $r_s$ (sound horizon) and $D_A$ (angular diameter distance) using the integration results.
    -   It computes the resulting $\Omega_m, \Omega_{DE}$, and $\theta_s$.
5.  **Residual Calculation**: It compares the calculated values to the targets and returns a "residual" (error).
6.  **Iteration**: The optimizer adjusts the shooting variables and repeats steps 2-5 until the error is minimized (near zero).

---

## 4. Verification of Physics and Code

### Physics Check
- **Crossing $w=-1$**: The combination of a standard field ($\phi$) and a phantom field ($\sigma$) allows the total dark energy equation of state $w_{DE}$ to effectively cross $-1$, which is the signature of Quintom models.
- **Background Evolution**: The notebook correctly implements the back-reaction of the scalar fields on the Hubble rate $H = \sqrt{\rho_{tot}}$.
- **Recombination**: The use of `calculate_z_dec_peebles` ensures the sound horizon is calculated at the correct physical epoch, rather than relying on a fixed $z=1089$ approximation.

### Code Check
- **Numerical Stability**: The ODE system uses the `Radau` method, which is suited for "stiff" equations (common in dark energy models where fields change scales rapidly).
- **Precision**: The integration uses high relative and absolute tolerances (`rtol=1e-12`, `atol=1e-14`), and `Simpson's` rule for integral quantities like $\theta_s$.
- **Unit Consistency**: The code works in "CLASS units" where $c=1$, and densities are scaled by $H_{0,class}^2$, ensuring compatibility with the main CLASS infrastructure.

### Conclusions
The implementation in `model_good2.ipynb` is physically sound and numerically robust. The shooting method successfully finds the unique set of $V_0$, $K_m$, and $h$ that satisfy the flat-universe constraints and CMB observations.
