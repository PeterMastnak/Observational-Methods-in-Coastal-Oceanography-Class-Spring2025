import numpy as np
import scipy.io as sio
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats

# Load the .mat file
data = sio.loadmat('data/X2DataFromIvy.mat')

# Extract Q (Dayflow) and X2 (CDECX2) data
Q = data['Dayflow'].flatten()  # Flow rate in m³/s
X2 = data['CDECX2'].flatten()  # Distance in km

# Define power law function X2 = AQ^n
def power_law(Q, A, n):
    return A * Q**n

# Fit power law
popt, pcov = curve_fit(power_law, Q, X2)
A, n = popt

# Calculate 95% confidence intervals
perr = np.sqrt(np.diag(pcov))
ci = 1.96 * perr  # 95% confidence interval

# Calculate standard error in X2
X2_fit = power_law(Q, A, n)
residuals = X2 - X2_fit
std_error = np.sqrt(np.mean(residuals**2))

# Perform logarithmic fit log10(X2) = a*log10(Q) + b
log_Q = np.log10(Q)
log_X2 = np.log10(X2)
slope, intercept, r_value, p_value, std_err = stats.linregress(log_Q, log_X2)
a = slope
b = intercept

# Print results
print(f"\nPower Law Fit Results (X2 = AQ^n):")
print(f"A = {A:.4f} ± {ci[0]:.4f} (95% CI)")
print(f"n = {n:.4f} ± {ci[1]:.4f} (95% CI)")
print(f"Standard Error in X2: {std_error:.4f} km")

print(f"\nLogarithmic Fit Results (log10(X2) = a*log10(Q) + b):")
print(f"a = {a:.4f}")
print(f"b = {b:.4f}")

print(f"\nComparison of exponents:")
print(f"Power law n: {n:.4f}")
print(f"Log fit a: {a:.4f}")
print(f"Difference: {abs(n-a):.4f}")

# Create plots
plt.figure(figsize=(12, 5))

# Plot 1: Direct relationship
plt.subplot(121)
plt.scatter(Q, X2, alpha=0.5, label='Data')
Q_smooth = np.linspace(min(Q), max(Q), 100)
plt.plot(Q_smooth, power_law(Q_smooth, A, n), 'r-', label='Power Law Fit')
plt.xlabel('Q (m³/s)')
plt.ylabel('X2 (km)')
plt.title('Power Law Fit: X2 = AQ^n')
plt.legend()
plt.grid(True)

# Plot 2: Log-log relationship
plt.subplot(122)
plt.scatter(log_Q, log_X2, alpha=0.5, label='Data')
plt.plot(log_Q, a*log_Q + b, 'r-', label='Linear Fit')
plt.xlabel('log10(Q)')
plt.ylabel('log10(X2)')
plt.title('Logarithmic Fit: log10(X2) = a*log10(Q) + b')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('Q_X2_analysis.png')
plt.close()

# Additional comments
print("\nAdditional Comments:")
print("1. The power law exponent (n) and logarithmic slope (a) should theoretically be the same")
print(f"2. The R-squared value for the log fit is: {r_value**2:.4f}")
print("3. The relationship appears to be well-described by a power law")
if abs(n-a) < 0.01:
    print("4. The two fitting methods give very consistent results")
else:
    print("4. There are some differences between the two fitting methods") 