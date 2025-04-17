import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def calculate_cd_finite_difference(xi2, xi0, x2, x0, h, u, g=9.81):
    """
    Calculate C_D using proper finite difference form
    
    Parameters:
    -----------
    xi2 : float
        Water surface elevation at position i+2 (m)
    xi0 : float
        Water surface elevation at position i (m)
    x2 : float
        Position of second measurement (m)
    x0 : float
        Position of first measurement (m)
    h : float
        Water depth at position i+1 (m)
    u : float
        Depth-averaged velocity at position i+1 (m/s)
    g : float
        Gravitational acceleration (m/s²), default 9.81
    
    Returns:
    --------
    cd : float
        Drag coefficient
    """
    # From finite difference equation: g(ξ_{i+2}-ξ_i)/(x_{i+2}-x_i) = -C_D*U_{i+1}|U_{i+1}|/h_{i+1}
    # Therefore: C_D = -g*h*(ξ_{i+2}-ξ_i)/(x_{i+2}-x_i)/(U_{i+1}|U_{i+1}|)
    return -g * h * (xi2 - xi0) / (x2 - x0) / (u * abs(u))

def calculate_cd_uncertainty(xi2, xi0, x2, x0, h, u, 
                           xi_uncertainty=0.001,  # 1mm uncertainty in surface elevation
                           x_uncertainty=0.01,    # 1cm uncertainty in position
                           h_uncertainty=0.01,    # 1cm uncertainty in depth
                           u_uncertainty=0.05):   # 5cm/s uncertainty in velocity
    """
    Calculate C_D and its uncertainty using finite difference approach
    
    Parameters:
    -----------
    xi2, xi0 : float
        Water surface elevations at positions i+2 and i (m)
    x2, x0 : float
        Positions of measurements (m)
    h : float
        Water depth at position i+1 (m)
    u : float
        Depth-averaged velocity at position i+1 (m/s)
    """
    g = 9.81
    
    # Calculate C_D using finite differences
    cd = calculate_cd_finite_difference(xi2, xi0, x2, x0, h, u)
    
    # Calculate partial derivatives for uncertainty propagation
    dx = x2 - x0
    dxi = xi2 - xi0
    
    # Partial derivatives with respect to each variable
    dcd_dxi2 = -g * h / (dx * u * abs(u))
    dcd_dxi0 = g * h / (dx * u * abs(u))
    dcd_dx2 = g * h * dxi / (dx**2 * u * abs(u))
    dcd_dx0 = -g * h * dxi / (dx**2 * u * abs(u))
    dcd_dh = -g * dxi / (dx * u * abs(u))
    dcd_du = 2 * g * h * dxi / (dx * u**2 * abs(u))
    
    # Calculate uncertainty using error propagation
    cd_uncertainty = np.sqrt(
        (dcd_dxi2 * xi_uncertainty)**2 +
        (dcd_dxi0 * xi_uncertainty)**2 +
        (dcd_dx2 * x_uncertainty)**2 +
        (dcd_dx0 * x_uncertainty)**2 +
        (dcd_dh * h_uncertainty)**2 +
        (dcd_du * u_uncertainty)**2
    )
    
    return cd, cd_uncertainty

def main():
    # Parameters
    h = 10  # water depth (m)
    g = 9.81  # gravitational acceleration (m/s²)
    
    # Create arrays for analysis
    separations = np.linspace(0.5, 5, 50)  # sensor separations from 0.5 to 5 meters
    velocities = np.linspace(0.1, 1.5, 50)  # velocities from 0.1 to 1.5 m/s
    
    # Initialize arrays for results
    cd_uncertainties = np.zeros((len(separations), len(velocities)))
    
    # Target C_D value for calculating expected surface elevation difference
    target_cd = 0.025  # middle of expected range (0.02-0.03)
    
    # Calculate uncertainties for each combination
    for i, sep in enumerate(separations):
        x0 = 0  # reference position
        x2 = sep  # separation distance
        
        for j, vel in enumerate(velocities):
            # Calculate expected surface elevation difference for target C_D
            # From rearranging finite difference equation:
            # (ξ_{i+2}-ξ_i) = (-C_D*U_{i+1}|U_{i+1}|*(x_{i+2}-x_i))/(g*h)
            xi0 = 0  # reference elevation
            xi2 = (-target_cd * vel * abs(vel) * sep) / (g * h)
            
            # Calculate uncertainty using finite differences
            cd, cd_uncertainty = calculate_cd_uncertainty(xi2, xi0, x2, x0, h, vel)
            cd_uncertainties[i, j] = cd_uncertainty
    
    # Create contour plot
    plt.figure(figsize=(12, 8))
    contour = plt.contourf(velocities, separations, cd_uncertainties, levels=20, cmap='viridis')
    plt.colorbar(contour, label='C_D Uncertainty')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Sensor Separation (m)')
    plt.title('Uncertainty in C_D using Finite Differences\nfor Different Velocities and Sensor Separations')
    
    # Add reference lines and annotations
    plt.axhline(y=2, color='r', linestyle='--', label='Recommended separation (2m)')
    plt.axvline(x=0.5, color='g', linestyle='--', label='Minimum reliable velocity (0.5 m/s)')
    plt.legend()
    
    plt.savefig('cd_uncertainty_contour.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create uncertainty vs separation plot for different velocities
    plt.figure(figsize=(12, 6))
    selected_velocities = [0.5, 1.0, 1.5]
    for vel in selected_velocities:
        idx = np.abs(velocities - vel).argmin()
        plt.plot(separations, cd_uncertainties[:, idx], label=f'v = {vel:.1f} m/s')
    
    plt.xlabel('Sensor Separation (m)')
    plt.ylabel('C_D Uncertainty')
    plt.title('C_D Uncertainty vs Sensor Separation using Finite Differences')
    plt.legend()
    plt.grid(True)
    plt.savefig('cd_uncertainty_vs_separation.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main() 