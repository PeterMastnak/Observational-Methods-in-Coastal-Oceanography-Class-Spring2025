#!/usr/bin/env python3
"""
Air-Sea Toolbox implementation for Monterey Bay
Calculate sensible and latent heat fluxes and shear stress
for a summer afternoon in Monterey Bay
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants
GRAVITY = 9.81  # m/s^2
CP_AIR = 1004.67  # J/kg/K, specific heat capacity of air at constant pressure
VON_KARMAN = 0.4  # von Karman constant
CP_SW = 3850.0  # J/kg/K, specific heat capacity of seawater
RHO_AIR = 1.22  # kg/m^3, density of air at sea level
RHO_SW = 1025.0  # kg/m^3, density of seawater
LE = 2.5e6  # J/kg, latent heat of evaporation

# --------- Air-Sea Bulk Flux Algorithms ---------

def qsat(T, P=1013.25):
    """
    Calculate saturation specific humidity (kg/kg)
    
    Parameters:
    ----------
    T : float or array_like
        Temperature in degrees Celsius
    P : float or array_like
        Pressure in mb (hPa)
        
    Returns:
    -------
    q : float or array_like
        Saturation specific humidity (kg/kg)
    """
    # Calculate saturation vapor pressure (mb) using Bolton's formula
    es = 6.112 * np.exp((17.67 * T) / (T + 243.5))
    # Convert to specific humidity
    q = 0.622 * es / (P - 0.378 * es)
    # Apply salinity correction for ocean surface (98% of freshwater value)
    return q

def satvap(T):
    """
    Calculate saturation vapor pressure (mb)
    
    Parameters:
    ----------
    T : float or array_like
        Temperature in degrees Celsius
        
    Returns:
    -------
    es : float or array_like
        Saturation vapor pressure (mb)
    """
    # Use Bolton's formula for saturation vapor pressure
    es = 6.112 * np.exp((17.67 * T) / (T + 243.5))
    return es

def rel_humidity_from_q(q, T, P=1013.25):
    """
    Calculate relative humidity from specific humidity
    
    Parameters:
    ----------
    q : float or array_like
        Specific humidity (kg/kg)
    T : float or array_like
        Temperature in degrees Celsius
    P : float or array_like
        Pressure in mb (hPa)
        
    Returns:
    -------
    rh : float or array_like
        Relative humidity (%)
    """
    qs = qsat(T, P)
    rh = 100.0 * q / qs
    return rh

def visc_air(Ta):
    """
    Calculate the kinematic viscosity of air
    
    Parameters:
    ----------
    Ta : float or array_like
        Air temperature in degrees Celsius
        
    Returns:
    -------
    visc : float or array_like
        Kinematic viscosity of air (m^2/s)
    """
    visc = 1.326e-5 * (1 + 6.542e-3 * Ta + 8.301e-6 * Ta**2 - 4.84e-9 * Ta**3)
    return visc

def calc_L(usr, tsr, qsr, Ta, theta_v):
    """
    Calculate the Monin-Obukhov length
    
    Parameters:
    ----------
    usr : float
        Friction velocity (m/s)
    tsr : float
        Temperature scaling parameter (K)
    qsr : float
        Humidity scaling parameter (kg/kg)
    Ta : float
        Air temperature (K)
    theta_v : float
        Virtual potential temperature (K)
        
    Returns:
    -------
    L : float
        Monin-Obukhov length (m)
    """
    # Virtual temperature scaling parameter
    tvsr = tsr * (1 + 0.61 * qsr)
    # Monin-Obukhov length
    L = (usr**2 * theta_v) / (VON_KARMAN * GRAVITY * tvsr)
    return L

def psim_stable(zeta):
    """
    Stability function for momentum in stable conditions
    
    Parameters:
    ----------
    zeta : float or array_like
        z/L, stability parameter
        
    Returns:
    -------
    psim : float or array_like
        Stability function for momentum
    """
    # Limit stability parameter to avoid divergence
    zeta = np.minimum(zeta, 10)
    psim = -5.0 * zeta
    return psim

def psim_unstable(zeta):
    """
    Stability function for momentum in unstable conditions
    
    Parameters:
    ----------
    zeta : float or array_like
        z/L, stability parameter
        
    Returns:
    -------
    psim : float or array_like
        Stability function for momentum
    """
    x = (1.0 - 16.0 * zeta)**0.25
    psim = 2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x*x) / 2.0) - 2.0 * np.arctan(x) + np.pi/2.0
    return psim

def psih_stable(zeta):
    """
    Stability function for heat in stable conditions
    
    Parameters:
    ----------
    zeta : float or array_like
        z/L, stability parameter
        
    Returns:
    -------
    psih : float or array_like
        Stability function for heat
    """
    # Limit stability parameter to avoid divergence
    zeta = np.minimum(zeta, 10)
    psih = -5.0 * zeta
    return psih

def psih_unstable(zeta):
    """
    Stability function for heat in unstable conditions
    
    Parameters:
    ----------
    zeta : float or array_like
        z/L, stability parameter
        
    Returns:
    -------
    psih : float or array_like
        Stability function for heat
    """
    x = (1.0 - 16.0 * zeta)**0.25
    psih = 2.0 * np.log((1.0 + x*x) / 2.0)
    return psih

def calculate_z0(usr, Ta):
    """
    Calculate roughness length for momentum using Charnock relation
    
    Parameters:
    ----------
    usr : float
        Friction velocity (m/s)
    Ta : float
        Air temperature in degrees Celsius
        
    Returns:
    -------
    z0 : float
        Roughness length for momentum (m)
    """
    # Charnock parameter
    alpha_c = 0.011
    # Smooth flow roughness length
    nu = visc_air(Ta)
    z0s = 0.11 * nu / usr
    # Roughness length from Charnock relation
    z0c = alpha_c * usr**2 / GRAVITY
    # Combined roughness length (Smith 1988)
    z0 = z0s + z0c
    return z0

def calculate_z0t(z0, Ta, usr):
    """
    Calculate roughness length for temperature
    
    Parameters:
    ----------
    z0 : float
        Roughness length for momentum (m)
    Ta : float
        Air temperature in degrees Celsius
    usr : float
        Friction velocity (m/s)
        
    Returns:
    -------
    z0t : float
        Roughness length for temperature (m)
    """
    # Reynolds number for roughness length
    nu = visc_air(Ta)
    Re = z0 * usr / nu
    # Roughness Reynolds number threshold
    Rec = 0.11
    
    if Re <= Rec:
        z0t = z0 * np.exp(0)  # Smooth flow
    else:
        z0t = z0 * np.exp(-2.67)  # Rough flow
    
    return z0t

def calculate_z0q(z0, Ta, usr):
    """
    Calculate roughness length for humidity
    
    Parameters:
    ----------
    z0 : float
        Roughness length for momentum (m)
    Ta : float
        Air temperature in degrees Celsius
    usr : float
        Friction velocity (m/s)
        
    Returns:
    -------
    z0q : float
        Roughness length for humidity (m)
    """
    # Reynolds number for roughness length
    nu = visc_air(Ta)
    Re = z0 * usr / nu
    # Roughness Reynolds number threshold
    Rec = 0.11
    
    if Re <= Rec:
        z0q = z0 * np.exp(0)  # Smooth flow
    else:
        z0q = z0 * np.exp(-2.67)  # Rough flow
    
    return z0q

def air_sea_fluxes(u, zu, SST, Ta, za, q, zq, Pa=1015.0, max_iter=20, min_usr=0.1):
    """
    Calculate air-sea fluxes using an iterative approach
    following COARE 3.0 algorithm
    
    Parameters:
    ----------
    u : float
        Wind speed (m/s) at height zu
    zu : float
        Height of wind measurement (m)
    SST : float
        Sea surface temperature (°C)
    Ta : float
        Air temperature (°C) at height za
    za : float
        Height of temperature measurement (m)
    q : float
        Specific humidity (kg/kg) at height zq
    zq : float
        Height of humidity measurement (m)
    Pa : float, optional
        Surface pressure (mb), default 1015.0
    max_iter : int, optional
        Maximum number of iterations, default 20
    min_usr : float, optional
        Minimum friction velocity (m/s), default 0.1
        
    Returns:
    -------
    dict
        Dictionary containing the flux components and other variables
    """
    # Convert temperatures to Kelvin for calculations
    Ta_K = Ta + 273.15
    SST_K = SST + 273.15
    
    # Calculate saturated specific humidity at sea surface
    qs = qsat(SST, Pa)
    
    # Initial neutral transfer coefficients
    Cdn_10 = 1.3e-3  # Neutral drag coefficient at 10m
    Chn_10 = 1.0e-3  # Neutral transfer coefficient for heat at 10m
    Cen_10 = 1.0e-3  # Neutral transfer coefficient for moisture at 10m
    
    # Initial guesses
    usr = max(min_usr, np.sqrt(Cdn_10 * u**2))  # Friction velocity
    tsr = 0.0     # Temperature scaling parameter
    qsr = 0.0     # Humidity scaling parameter
    z0 = 0.0001  # Roughness length initial guess
    z0t = 0.0001 # Thermal roughness length initial guess
    z0q = 0.0001 # Humidity roughness length initial guess
    
    # Calculate virtual temperature
    theta_v = Ta_K * (1.0 + 0.61 * q)
    
    # Monin-Obukhov length initial guess (neutral)
    L = 1000.0  # Large positive value for near-neutral conditions
    
    # Iterative calculation
    for i in range(max_iter):
        # Stability parameter
        zeta_u = zu / L
        zeta_a = za / L
        zeta_q = zq / L
        
        # Stability functions
        if L > 0:  # Stable
            psim_u = psim_stable(zeta_u)
            psih_a = psih_stable(zeta_a)
            psih_q = psih_stable(zeta_q)
        else:      # Unstable
            psim_u = psim_unstable(zeta_u)
            psih_a = psih_unstable(zeta_a)
            psih_q = psih_unstable(zeta_q)
        
        # Update fluxes
        usr = u * VON_KARMAN / (np.log(zu/z0) - psim_u)
        usr = max(usr, min_usr)  # Apply minimum value to avoid division by zero
        
        # Update roughness lengths based on wind stress
        z0 = calculate_z0(usr, Ta)
        z0t = calculate_z0t(z0, Ta, usr)
        z0q = calculate_z0q(z0, Ta, usr)
        
        # Calculate transfer coefficients
        Cd = (VON_KARMAN / (np.log(zu/z0) - psim_u))**2
        Ch = (VON_KARMAN / (np.log(zu/z0) - psim_u)) * (VON_KARMAN / (np.log(za/z0t) - psih_a))
        Ce = (VON_KARMAN / (np.log(zu/z0) - psim_u)) * (VON_KARMAN / (np.log(zq/z0q) - psih_q))
        
        # Update scaling parameters (positive for upward flux, from ocean to atmosphere)
        tsr = Ch / Cd * (Ta - SST)  # Changed sign to match convention
        qsr = Ce / Cd * (q - qs)    # Changed sign to match convention
        
        # Update Monin-Obukhov length
        L_old = L
        L = calc_L(usr, tsr, qsr, Ta_K, theta_v)
        
        # Check for convergence
        if abs(L - L_old) < 0.001 * abs(L_old):
            break
    
    # Calculate fluxes (positive upward, from ocean to atmosphere)
    tau = RHO_AIR * usr**2  # Stress (N/m^2)
    Hs = -RHO_AIR * CP_AIR * usr * tsr  # Sensible heat flux (W/m^2)
    Hl = -RHO_AIR * LE * usr * qsr  # Latent heat flux (W/m^2)
    
    # Return results
    return {
        'tau': tau,          # Wind stress (N/m^2)
        'Hs': Hs,            # Sensible heat flux (W/m^2)
        'Hl': Hl,            # Latent heat flux (W/m^2)
        'usr': usr,          # Friction velocity (m/s)
        'tsr': tsr,          # Temperature scaling parameter (K)
        'qsr': qsr,          # Humidity scaling parameter (kg/kg)
        'z0': z0,            # Roughness length for momentum (m)
        'z0t': z0t,          # Roughness length for temperature (m)
        'z0q': z0q,          # Roughness length for humidity (m)
        'L': L,              # Monin-Obukhov length (m)
        'Cd': Cd,            # Drag coefficient
        'Ch': Ch,            # Transfer coefficient for heat
        'Ce': Ce             # Transfer coefficient for moisture
    }

# --------- Main Program ---------

def main():
    """
    Main function to calculate and display air-sea fluxes for Monterey Bay
    in summer afternoon conditions
    """
    print("Air-Sea Flux Calculations for Monterey Bay (Summer Afternoon)")
    print("-" * 60)
    
    # Input values for a summer afternoon in Monterey Bay
    # Based on NOAA buoy data for Monterey Bay in summer
    u = 5.0        # Wind speed (m/s) - moderate afternoon sea breeze
    zu = 2.0       # Height of wind measurement (m) above water surface
    SST = 15.5     # Sea surface temperature (°C) - typical summer value
    Ta = 17.0      # Air temperature (°C) - warmer than sea surface due to land warming
    za = 2.0       # Height of temperature measurement (m)
    RH = 80.0      # Relative humidity (%) - typical marine environment
    zq = 2.0       # Height of humidity measurement (m)
    Pa = 1015.0    # Atmospheric pressure (mb)
    
    # Calculate specific humidity from relative humidity
    q_sat_air = qsat(Ta, Pa)
    q = RH / 100.0 * q_sat_air
    
    # Calculate saturated specific humidity at sea surface
    qs = qsat(SST, Pa)
    
    # Print humidity gradient information for debugging
    dq = q - qs
    print(f"Input Parameters:")
    print(f"  Wind speed at {zu}m: {u:.1f} m/s")
    print(f"  Sea surface temperature: {SST:.1f} °C")
    print(f"  Air temperature at {za}m: {Ta:.1f} °C")
    print(f"  Relative humidity at {zq}m: {RH:.1f} %")
    print(f"  Specific humidity at {zq}m: {q*1000:.2f} g/kg")
    print(f"  Saturated specific humidity at sea surface: {qs*1000:.2f} g/kg")
    print(f"  Specific humidity difference (air-sea): {dq*1000:.2f} g/kg")
    print(f"  Atmospheric pressure: {Pa:.1f} mb")
    print()
    
    # Calculate air-sea fluxes
    fluxes = air_sea_fluxes(u, zu, SST, Ta, za, q, zq, Pa)
    
    # IMPORTANT: Sign convention for fluxes
    # In our implementation:
    # - Positive heat flux = heat transfer from ocean to atmosphere (ocean cooling)
    # - Negative heat flux = heat transfer from atmosphere to ocean (ocean warming)
    # This is the oceanographic convention (positive upward)
    Hs = fluxes['Hs']
    Hl = fluxes['Hl']
    
    # Display results
    print("Calculated Fluxes:")
    print(f"  Wind stress (τ): {fluxes['tau']*1000:.2f} mN/m²")
    print(f"  Sensible heat flux (Hs): {Hs:.2f} W/m²")
    print(f"  Latent heat flux (Hl): {Hl:.2f} W/m²")
    print(f"  Total turbulent heat flux (Hs + Hl): {Hs + Hl:.2f} W/m²")
    
    # Interpretation of fluxes
    if Hs > 0:
        print(f"  -> Sensible heat is transferring from ocean to atmosphere (ocean cooling)")
    else:
        print(f"  -> Sensible heat is transferring from atmosphere to ocean (ocean warming)")
        
    if Hl > 0:
        print(f"  -> Latent heat is transferring from ocean to atmosphere (evaporation)")
    else:
        print(f"  -> Latent heat is transferring from atmosphere to ocean (condensation)")
    print()
    
    print("Transfer Coefficients:")
    print(f"  Drag coefficient (Cd): {fluxes['Cd']*1000:.2f} × 10⁻³")
    print(f"  Heat transfer coefficient (Ch): {fluxes['Ch']*1000:.2f} × 10⁻³")
    print(f"  Moisture transfer coefficient (Ce): {fluxes['Ce']*1000:.2f} × 10⁻³")
    print()
    
    print("Other Parameters:")
    print(f"  Friction velocity (u*): {fluxes['usr']:.3f} m/s")
    print(f"  Roughness length (z₀): {fluxes['z0']*1000:.3f} mm")
    print(f"  Monin-Obukhov length (L): {fluxes['L']:.1f} m")
    stability = "stable" if fluxes['L'] > 0 else "unstable"
    print(f"  Atmospheric stability: {stability}")
    
    # Calculate evaporation rate (correct sign if needed)
    latent_heat_flux = abs(Hl) if Hl > 0 else 0  # Use only positive latent heat flux for evaporation
    evap_rate = latent_heat_flux / LE  # kg/m²/s
    evap_mm_per_day = evap_rate * 86400  # mm/day
    
    print(f"  Evaporation rate: {evap_mm_per_day:.2f} mm/day")
    
    # Create simple plot of results
    plot_results(fluxes)

def plot_results(fluxes):
    """Create simple plots of the results"""
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    # Heat fluxes plot
    fluxes_names = ['Sensible Heat', 'Latent Heat']
    fluxes_values = [fluxes['Hs'], fluxes['Hl']]
    
    # Add total as third bar
    fluxes_names.append('Total Heat')
    fluxes_values.append(fluxes['Hs'] + fluxes['Hl'])
    
    # Create bar chart with conditional colors (positive/negative)
    colors = ['blue' if x < 0 else 'red' for x in fluxes_values]
    axs[0].bar(fluxes_names, fluxes_values, color=colors, alpha=0.7)
    axs[0].set_ylabel('Heat Flux (W/m²)')
    axs[0].set_title('Turbulent Heat Fluxes\n(positive = ocean to atmosphere)')
    axs[0].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Add flux values as text
    for i, v in enumerate(fluxes_values):
        axs[0].text(i, v + np.sign(v)*5, f"{v:.1f}", ha='center')
    
    # Transfer coefficients plot
    coef_names = ['Drag (Cd)', 'Heat (Ch)', 'Moisture (Ce)']
    coef_values = [fluxes['Cd']*1000, fluxes['Ch']*1000, fluxes['Ce']*1000]  # Convert to 10^-3
    
    axs[1].bar(coef_names, coef_values, color='green', alpha=0.7)
    axs[1].set_ylabel('Coefficient (×10⁻³)')
    axs[1].set_title('Transfer Coefficients')
    
    # Add coefficient values as text
    for i, v in enumerate(coef_values):
        axs[1].text(i, v + 0.05, f"{v:.2f}", ha='center')
    
    plt.tight_layout()
    plt.savefig('monterey_bay_fluxes.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as 'monterey_bay_fluxes.png'")
    plt.show()

if __name__ == "__main__":
    main() 