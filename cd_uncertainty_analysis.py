import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def calculate_cd_uncertainty(delta_p, rho, h, u, delta_h_uncertainty=0.01, delta_p_uncertainty=0.001):
    """
    Calculate C_D and its uncertainty based on pressure difference measurements
    
    Parameters:
    -----------
    delta_p : float
        Pressure difference between sensors (Pa)
    rho : float
        Water density (kg/m³)
    h : float
        Water depth (m)
    u : float
        Depth-averaged velocity (m/s)
    delta_h_uncertainty : float
        Uncertainty in depth measurement (m)
    delta_p_uncertainty : float
        Uncertainty in pressure difference measurement (Pa)
    
    Returns:
    --------
    cd : float
        Drag coefficient
    cd_uncertainty : float
        Uncertainty in drag coefficient
    """
    g = 9.81  # gravitational acceleration (m/s²)
    
    # Calculate C_D
    cd = 2 * delta_p / (rho * u**2)
    
    # Calculate partial derivatives for uncertainty propagation
    dcd_dp = 2 / (rho * u**2)
    dcd_du = -4 * delta_p / (rho * u**3)
    dcd_dh = 0  # C_D is not directly dependent on h
    
    # Calculate uncertainty
    cd_uncertainty = np.sqrt(
        (dcd_dp * delta_p_uncertainty)**2 +
        (dcd_du * delta_h_uncertainty)**2 +
        (dcd_dh * delta_h_uncertainty)**2
    )
    
    return cd, cd_uncertainty

def main():
    # Parameters
    rho = 1025  # seawater density (kg/m³)
    h = 10  # water depth (m)
    
    # Create arrays for analysis
    separations = np.linspace(0.5, 5, 50)  # sensor separations from 0.5 to 5 meters
    velocities = np.linspace(0.1, 1.5, 50)  # velocities from 0.1 to 1.5 m/s
    
    # Initialize arrays for results
    cd_uncertainties = np.zeros((len(separations), len(velocities)))
    
    # Calculate uncertainties for each combination
    for i, sep in enumerate(separations):
        for j, vel in enumerate(velocities):
            # Calculate pressure difference (simplified)
            delta_p = 0.5 * rho * vel**2 * sep / h
            cd, cd_uncertainty = calculate_cd_uncertainty(delta_p, rho, h, vel)
            cd_uncertainties[i, j] = cd_uncertainty
    
    # Create contour plot
    plt.figure(figsize=(10, 8))
    plt.contourf(velocities, separations, cd_uncertainties, levels=20)
    plt.colorbar(label='C_D Uncertainty')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Sensor Separation (m)')
    plt.title('Uncertainty in C_D for Different Velocities and Sensor Separations')
    
    # Add reference lines for typical C_D values
    plt.axhline(y=2, color='r', linestyle='--', label='Typical separation (2m)')
    
    plt.savefig('cd_uncertainty_contour.png')
    plt.close()
    
    # Create a plot showing uncertainty vs separation for different velocities
    plt.figure(figsize=(10, 6))
    selected_velocities = [0.5, 1.0, 1.5]
    for vel in selected_velocities:
        idx = np.abs(velocities - vel).argmin()
        plt.plot(separations, cd_uncertainties[:, idx], label=f'v = {vel:.1f} m/s')
    
    plt.xlabel('Sensor Separation (m)')
    plt.ylabel('C_D Uncertainty')
    plt.title('C_D Uncertainty vs Sensor Separation for Different Velocities')
    plt.legend()
    plt.grid(True)
    plt.savefig('cd_uncertainty_vs_separation.png')
    plt.close()

if __name__ == "__main__":
    main() 