from classy import Class
import numpy as np

def run_test(delta=0.0, coupled=False):
    cosmo = Class()
    params = {
        'H0': 67.55,
        'omega_b': 0.022032,
        'omega_cdm': 0.12038,
        'Omega_k': 0,
        'Omega_Lambda': 0,
        'Omega_fld': 0,
        'Omega_qtm': -1,
        'lambda_qtm': 1.0,
        'delta_qtm': delta,
        'coupled_baryon_qtm': coupled,
        'coupled_cdm_qtm': coupled,
        'input_verbose': 1,
        'background_verbose': 1
    }
    
    print(f"\n--- Running with delta_qtm={delta}, coupled={coupled} ---")
    try:
        cosmo.set(params)
        cosmo.compute(["background"])
        z_eq = cosmo.z_eq()
        print(f"CLASS z_eq: {z_eq}")
        
        # Check background data
        bg = cosmo.get_background()
        idx_mid = len(bg['z']) // 2
        z_val = bg['z'][idx_mid]
        rho_b = bg['(.)rho_b'][idx_mid]
        a = 1.0 / (1.0 + z_val)
        
        # Predicted rho_b if uncoupled (rho ~ a^-3)
        # We need the value at z=0 to compare correctly, but let's just check the scaling ratio if possible.
        # Or just use the first and last points.
        rho_b_0 = bg['(.)rho_b'][0]
        rho_b_expected_uncoupled = rho_b_0 / (a**3)
        
        ratio = rho_b / rho_b_expected_uncoupled
        print(f"At z={z_val:.2f}, rho_b / rho_b_standard = {ratio:.6f}")
        
        return z_eq, ratio
        
    except Exception as e:
        print(f"Error: {e}")
        return None, None
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()

# Reference LCDM-like (delta=0, coupled='no')
z_ref, ratio_ref = run_test(delta=0.0, coupled='no')

# Case with delta=0.1 but coupled='no'
z_nosupply, ratio_nosupply = run_test(delta=0.1, coupled='no')

# Case with delta=0.1 and coupled=True (Testing boolean fix)
z_coupled, ratio_coupled = run_test(delta=0.1, coupled=True)

print("\nSUMMARY:")
print(f"Reference z_eq: {z_ref}")
print(f"Delta 0.1 (No coupling) z_eq: {z_nosupply}")
print(f"Delta 0.1 (Coupled) z_eq: {z_coupled}")

if z_nosupply is not None and z_coupled is not None:
    if abs(z_nosupply - z_ref) < 1e-3:
        print("SUCCESS: delta_qtm=0.1 with coupled='no' correctly results in LCDM-like z_eq.")
    else:
        print("FAILURE: delta_qtm=0.1 with coupled='no' shifted z_eq unexpectedly.")
        
    if abs(z_coupled - z_ref) > 1.0:
        print("SUCCESS: delta_qtm=0.1 with coupled='yes' correctly shifted z_eq.")
    else:
        print("FAILURE: delta_qtm=0.1 with coupled='yes' did NOT shift z_eq significantly.")
