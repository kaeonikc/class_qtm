import numpy as np
from scipy.integrate import solve_ivp
import scipy.constants as const
from numba import jit

# --- High Precision Constants ---
# CODATA 2018 / SI
m_e = const.m_e        # Electron mass [kg]
k_B = const.k          # Boltzmann constant [J/K]
hbar = const.hbar      # Reduced Planck constant [J s]
c = const.c            # Speed of light in vacuum [m/s]
sigma_T = const.physical_constants['Thomson cross section'][0] # [m^2]
m_p = const.m_p        # Proton mass [kg]
eV = const.eV          # Electron volt [J]
sigma_SB = const.sigma # Stefan-Boltzmann

# Energy levels
E_I = 13.598434005136 * eV  # Ionization energy of H (13.6 eV) [J]
E_Lya = 10.19881 * eV       # Lyman-alpha energy (~10.2 eV) [J]
Lambda_2s = 8.22458         # 2s -> 1s two-photon decay rate [s^-1]

# Planck constant for Numba (const.h is a float, so Numba can capture it)
h_planck = const.h

# --- JIT Compiled ODE Physics ---
# We pass all constants as arguments to avoid closure overhead and allow Numba compilation
@jit(nopython=True, cache=True)
def calc_derivs(z, Xe, Tm, 
                n_H0, Omega_m0, Omega_rad0, H0_SI, Tcmb0, F_factor):
    
    # Clamp physical bounds
    if Xe > 1.0: Xe = 1.0
    if Xe < 1e-12: Xe = 1e-12
    
    one_plus_z = 1.0 + z
    Tr = Tcmb0 * one_plus_z
    nH = n_H0 * one_plus_z**3
    
    # Hubble
    E_z = np.sqrt(Omega_m0 * one_plus_z**3 + Omega_rad0 * one_plus_z**4)
    H_z = H0_SI * E_z
    
    # Matter Temp Evolution
    comp_coeff = (32.0 * sigma_T * sigma_SB * Tr**4) / (3.0 * m_e * c**2)
    dTm_dt = -2.0 * H_z * Tm + comp_coeff * (Tr - Tm) * Xe
    
    # dt/dz = -1 / (H_z * (1+z))
    driver = -1.0 / (H_z * one_plus_z)
    dTm_dz = dTm_dt * driver
    
    # Recombination Rate
    t = Tm / 1e4
    # Standard Recfast alpha is F * alpha. 
    # But for Peebles C_r regime, we apply F to C_r to simulate multi-level speedup.
    a_c = 4.309; b_c = -0.6166; c_c = 0.6703; d_c = 0.5300
    alpha_B = 1e-19 * (a_c * t**b_c) / (1.0 + c_c * t**d_c)
    
    # Beta calculation (consistent with alpha)
    # k_B, m_e, h_planck etc are globals, Numba captures them as constants if they are simple scalars.
    c_factor = ((2.0 * np.pi * m_e * k_B * Tm) / (h_planck**2))**1.5
    beta_2 = alpha_B * c_factor * np.exp(-(E_I/4.0) / (k_B * Tm))
    beta_ground = beta_2 * np.exp(-E_Lya / (k_B * Tm))
    
    lambda_Lya = (2.0 * np.pi * hbar * c) / E_Lya 
    
    if Xe >= 1.0:
        C_r = 1.0 
    else:
        K_factor = (nH * (1.0-Xe) * lambda_Lya**3) / (8.0 * np.pi * H_z)
        Lambda_alpha = 1.0 / K_factor
        C_r = (Lambda_alpha + Lambda_2s) / (Lambda_alpha + Lambda_2s + beta_2)
    
    # Apply Fudge Factor to C_r (Effective escape probability enhancer)
    # Calibrated: F ~ 1.1 increases z_dec to ~1090
    C_r *= F_factor
    
    # RHS
    RHS_t = - C_r * (alpha_B * nH * Xe**2 - beta_ground * (1.0 - Xe))
    dXe_dz = RHS_t * driver
    
    return dXe_dz, dTm_dz

def calculate_z_dec_peebles(cosmic_params, return_tables=False):
    """
    Calculate the decoupling redshift using Peebles recombination model (Recfast-style physics).
    Optimized with Numba for ODE and Numpy for vectorized Saha calculation.

    This function calculates two definitions of recombination redshift:
    1. z_dec (Peak): The peak of the visibility function g(z). Typically around 1077.
    2. z_star (Tau=1): The redshift where optical depth to Thomson scattering is unity. Typically around 1090.

    Parameters:
    -----------
    cosmic_params : dict
        Dictionary containing cosmological parameters:
        - 'h': dimensionless Hubble parameter (H0 / 100 km/s/Mpc)
        - 'Omega_b0': Baryon density parameter (Omega_b = rho_b / rho_crit). approx 0.05.
                      (NOT physical density omega_b = Omega_b h^2 approx 0.022).
        - 'Omega_m0': Matter density parameter (Omega_m). approx 0.3.
        - 'Omega_rad0': Radiation density parameter.
        - 'Tcmb0': CMB temperature today [K].
        - 'Yp': Helium mass fraction.
        - 'N_eff': Effective neutrino species.
        - 'fudge_factor': Peebles coefficient multiplier.
    return_tables : bool
        If True, returns a dictionary with full evolution tables.
    """

    # Extract parameters
    h = cosmic_params['h']
    Omega_b0 = cosmic_params['Omega_b0']
    Omega_m0 = cosmic_params['Omega_m0']
    Omega_rad0 = cosmic_params.get('Omega_rad0', 0.0)
    Tcmb0 = cosmic_params.get('Tcmb0', 2.7255)
    Yp = cosmic_params.get('Yp', 0.245)
    N_eff = cosmic_params.get('N_eff', 3.046)
    
    # Fudge Factor for Peebles coefficient
    # Calibrated to F=1.180 (Numba/LSODA) to yield z_dec ~ 1090
    F_factor = cosmic_params.get('fudge_factor', 1.140)

    # Derived Constants
    Mpc_in_m = const.parsec * 1e6 
    H0_val_kmsMpc = 100.0 * h
    H0_SI = (H0_val_kmsMpc * 1000.0) / Mpc_in_m 
    
    # Calculate Omega_rad0 consistency if not provided
    if Omega_rad0 is None or Omega_rad0 == 0.0:
        rho_crit0_SI = 3.0 * H0_SI**2 / (8.0 * np.pi * const.G) 
        rho_gamma = 4.0 * const.sigma * Tcmb0**4 / const.c**3 
        nu_factor = 1.0 + (7.0/8.0)*(4.0/11.0)**(4.0/3.0) * N_eff
        rho_rad_val = rho_gamma * nu_factor
        Omega_rad0 = rho_rad_val / rho_crit0_SI
    
    rho_crit0_SI = 3.0 * H0_SI**2 / (8.0 * np.pi * const.G)
    m_H = m_p + m_e
    n_H0 = ((1.0 - Yp) * Omega_b0 * rho_crit0_SI) / m_H 

    # --- Wrapped ODE for Scipy ---
    # Scipy expects f(t, y)
    def ode_wrapper(z, state):
        dXe, dTm = calc_derivs(z, state[0], state[1], 
                               n_H0, Omega_m0, Omega_rad0, H0_SI, Tcmb0, F_factor)
        return [dXe, dTm]

    # --- Vectorized Saha Check ---
    # Replaces the 10,000 loop
    z_start = 1600.0
    z_end = 10.0
    z_scan = np.linspace(z_start, z_end, 10000) # Descending array
    
    Tr_arr = Tcmb0 * (1.0 + z_scan)
    Tm_arr = Tr_arr # Saha approx
    nH_arr = n_H0 * (1.0 + z_scan)**3
    
    t_arr = Tm_arr / 1e4
    a_c = 4.309; b_c = -0.6166; c_c = 0.6703; d_c = 0.5300
    alpha_B_arr = 1e-19 * (a_c * t_arr**b_c) / (1.0 + c_c * t_arr**d_c)
    
    const_fac_saha = ((2.0 * np.pi * m_e * k_B * Tm_arr) / h_planck**2)**1.5
    beta_2_arr = alpha_B_arr * const_fac_saha * np.exp(-(E_I/4.0) / (k_B * Tm_arr))
    beta_ground_arr = beta_2_arr * np.exp(-E_Lya / (k_B * Tm_arr))
    
    S_saha = beta_ground_arr / (alpha_B_arr * nH_arr)
    Xe_saha = (-S_saha + np.sqrt(S_saha**2 + 4.0 * S_saha)) / 2.0
    Xe_saha = np.minimum(Xe_saha, 1.0)
    
    # Find transition point Xe < 0.99
    # We want the first index (highest z) where Xe drops below 0.99
    mask = Xe_saha < 0.99
    if np.any(mask):
        idx_switch = np.argmax(mask) # returns first True index
    else:
        idx_switch = len(z_scan) - 1 # Always Saha (unlikely), integrate ODE from start
        
    # Prepare result arrays for Saha regime
    z_saha_regime = z_scan[:idx_switch]
    Xe_saha_regime = Xe_saha[:idx_switch]
    Tm_saha_regime = Tm_arr[:idx_switch]
    
    # Initial conditions for ODE
    # The ODE should start from the last Saha point, or z_start if no Saha points
    if idx_switch > 0:
        z_ode_start = z_scan[idx_switch-1]
        Xe_ini = Xe_saha[idx_switch-1] 
        Tm_ini = Tm_arr[idx_switch-1]
    else: # If Xe_saha is already < 0.99 at z_start
        z_ode_start = z_scan[0]
        Xe_ini = Xe_saha[0]
        Tm_ini = Tm_arr[0]
    
    # Solve ODE
    sol = solve_ivp(ode_wrapper, (z_ode_start, z_end), [Xe_ini, Tm_ini], 
                    method='LSODA', rtol=1e-6, atol=1e-10)
    
    if not sol.success:
        print(f"ODE Warning: {sol.message}")
        
    # Merge Results
    z_ode_out = sol.t
    Xe_ode_out = sol.y[0]
    Tm_ode_out = sol.y[1]
    
    # Interpolate ODE results onto the z_scan global grid for z < z_switch
    # z_scan is descending. solve_ivp.t is also descending.
    z_ode_target = z_scan[idx_switch:]
    
    # Ensure interpolation arrays are sorted ascending for np.interp
    # z_ode_out is descending, z_ode_target is descending.
    # np.interp expects x to be increasing.
    Xe_ode_interp = np.interp(z_ode_target[::-1], z_ode_out[::-1], Xe_ode_out[::-1])[::-1]
    Tm_ode_interp = np.interp(z_ode_target[::-1], z_ode_out[::-1], Tm_ode_out[::-1])[::-1]
    
    z_grid = np.concatenate([z_saha_regime, z_ode_target])
    Xe_grid = np.concatenate([Xe_saha_regime, Xe_ode_interp])
    Tm_grid = np.concatenate([Tm_saha_regime, Tm_ode_interp])
    
    # --- Optical Depth & g(z) ---
    # Vectorized calculation on z_grid
    
    one_plus_z = 1.0 + z_grid
    nH_grid = n_H0 * one_plus_z**3
    ne_grid = Xe_grid * nH_grid
    Ez_grid = np.sqrt(Omega_m0 * one_plus_z**3 + Omega_rad0 * one_plus_z**4)
    Hz_grid = H0_SI * Ez_grid
    
    dtau_dz_grid = ne_grid * sigma_T * c / (one_plus_z * Hz_grid)
    
    # Optical Depth tau(z) = Integral_0^z (n_e sigma_T c / (H(z)(1+z))) dz
    # Our grid is z_grid (descending from z_start to z_end).
    # To integrate from z=0 up to z, we need to integrate dtau/dz from z_end up to z_start.
    # cumulative_trapezoid integrates from left to right.
    # So, reverse arrays, integrate, then reverse result.
    
    from scipy.integrate import cumulative_trapezoid
    
    # Flip to increasing z for integration
    z_inc = z_grid[::-1]
    dtau_inc = dtau_dz_grid[::-1]
    
    # Integrate from z_end (approx 0) up to z_start
    tau_inc = cumulative_trapezoid(dtau_inc, z_inc, initial=0)
    
    # Flip back to match original z_grid order (descending)
    tau_arr = tau_inc[::-1]
    
    # Visibility Function g(z)
    g_arr = dtau_dz_grid * np.exp(-tau_arr)
    
    # Find z_dec (Peak of g)
    idx_max = np.argmax(g_arr)
    z_dec_peak = z_grid[idx_max]
    
    # Find z_star (tau = 1.0)
    idx_tau1 = np.abs(tau_arr - 1.0).argmin()
    z_star = z_grid[idx_tau1]

    # --- Baryon Drag Epoch ---
    # R = 3 * rho_b / (4 * rho_gamma)
    # rho_b ~ (1+z)^3, rho_gamma ~ (1+z)^4
    # R = 3 * Omega_b0 / (4 * Omega_gamma0 * (1+z))
    # Note: Omega_rad0 includes neutrinos usually. We need Omega_gamma0 specifically.
    
    rho_crit0_SI = 3.0 * H0_SI**2 / (8.0 * np.pi * const.G) 
    rho_gamma0 = 4.0 * const.sigma * Tcmb0**4 / const.c**3 
    Omega_gamma0 = rho_gamma0 / rho_crit0_SI
    
    R_grid = 3.0 * Omega_b0 / (4.0 * Omega_gamma0 * (1.0 + z_grid))
    
    # dtau_drag/dz = dtau/dz / R
    dtau_drag_dz_grid = dtau_dz_grid / R_grid
    
    # Integrate to find tau_drag
    # Reverse to integrate from z=0 (approx) up to z_start
    dtau_drag_inc = dtau_drag_dz_grid[::-1]
    tau_drag_inc = cumulative_trapezoid(dtau_drag_inc, z_inc, initial=0)
    tau_drag_arr = tau_drag_inc[::-1]
    
    # Find z_drag (tau_drag = 1.0)
    idx_tau_drag1 = np.abs(tau_drag_arr - 1.0).argmin()
    z_drag = z_grid[idx_tau_drag1]
    
    if return_tables:
        return {
            "z_dec": z_dec_peak,
            "z_star": z_star,
            "z_grid": z_grid,
            "Xe_grid": Xe_grid,
            "Tm_grid": Tm_grid,
            "tau_grid": tau_arr,
            "g_grid": g_arr,
            "dtau_dz_grid": dtau_dz_grid,
            "R_grid": R_grid,
            "tau_drag_grid": tau_drag_arr,
            "z_drag": z_drag
        }
    else:
        return z_dec_peak

def calculate_z_drag(cosmic_params):
    """
    Convenience function to calculate the baryon drag epoch redshift.
    """
    res = calculate_z_dec_peebles(cosmic_params, return_tables=True)
    return res['z_drag']

if __name__ == "__main__":
    # Test Block
    params = {
        'h': 0.683,
        'Omega_b0': 0.02237 / (0.683**2),
        'Omega_m0': 0.31,
        'Omega_rad0': 9.1e-5,
        'Tcmb0': 2.7255
    }
    
    print("Running Peebles Recombination Test...")
    ombh2 = params['Omega_b0'] * params['h']**2
    ommh2 = params['Omega_m0'] * params['h']**2
    print(f"Inputs: h={params['h']}, Omega_b0={params['Omega_b0']:.4f}, Omega_m0={params['Omega_m0']:.4f}")
    print(f"Physical: omega_b={ombh2:.5f}, omega_m={ommh2:.5f}")
    
    res = calculate_z_dec_peebles(params, return_tables=True)
    
    print("-" * 40)
    print(f"z_dec (Peak Visibility):     {res['z_dec']:.4f}")
    print(f"z_star (Optical Depth=1):    {res['z_star']:.4f}")
    print(f"z_drag (Drag Depth=1):       {res['z_drag']:.4f}")
    print("-" * 40)
    print("Test Complete.")
