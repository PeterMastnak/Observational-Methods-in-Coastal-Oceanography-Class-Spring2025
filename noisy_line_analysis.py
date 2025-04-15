import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

def generate_noisy_data(x_noise_factor, y_noise_factor):
    # Generate base X data
    X = np.linspace(-1, 1, 201)  # Same as -1:0.01:1 in MATLAB
    nX = len(X)
    
    # Generate noise
    r1 = x_noise_factor * np.random.randn(nX)
    r2 = y_noise_factor * np.random.randn(nX)
    
    # Create noisy X and Y
    X_noisy = X + r1
    Y_noisy = X + r2  # True relationship is Y = X
    
    return X_noisy, Y_noisy

def fit_and_analyze(X, Y, title):
    # Fit Y = AX + B
    A, B = np.polyfit(X, Y, 1)
    y_pred = A * X + B
    
    # Calculate 95% confidence intervals
    n = len(X)
    x_mean = np.mean(X)
    sum_sq_x = np.sum((X - x_mean) ** 2)
    std_err = np.sqrt(np.sum((Y - y_pred) ** 2) / (n - 2))
    
    # Standard error of slope
    se_A = std_err / np.sqrt(sum_sq_x)
    ci_A = 1.96 * se_A  # 95% confidence interval
    
    # Fit X = A'Y + B' (reverse fit)
    A_prime, B_prime = np.polyfit(Y, X, 1)
    
    # Calculate R-squared
    r_squared = stats.pearsonr(X, Y)[0]**2
    
    print(f"\n{title}")
    print(f"Forward fit (Y = AX + B):")
    print(f"A = {A:.4f} ± {ci_A:.4f} (95% CI)")
    print(f"B = {B:.4f}")
    print(f"Reverse fit (X = A'Y + B'):")
    print(f"A' = {A_prime:.4f}")
    print(f"B' = {B_prime:.4f}")
    print(f"1/A = {1/A:.4f}")
    print(f"R² = {r_squared:.4f}")
    
    return A, B, A_prime, B_prime, r_squared

# Test different noise factor combinations
noise_combinations = [
    (0.0, 0.2, "Case 1: No X noise, moderate Y noise"),
    (0.2, 0.2, "Case 2: Equal X and Y noise"),
    (0.4, 0.2, "Case 3: More X noise than Y noise")
]

plt.figure(figsize=(15, 5))

for idx, (x_noise, y_noise, title) in enumerate(noise_combinations, 1):
    X, Y = generate_noisy_data(x_noise, y_noise)
    
    plt.subplot(1, 3, idx)
    plt.scatter(X, Y, alpha=0.5, label='Data')
    
    # Perform fits and get results
    A, B, A_prime, B_prime, r_squared = fit_and_analyze(X, Y, title)
    
    # Plot fits
    x_line = np.linspace(min(X), max(X), 100)
    plt.plot(x_line, A*x_line + B, 'r-', label=f'Y = AX + B')
    plt.plot(x_line, (x_line - B_prime)/A_prime, 'g--', label=f'X = A\'Y + B\'')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f"{title}\nR² = {r_squared:.4f}")
    plt.legend()
    plt.grid(True)

plt.tight_layout()
plt.savefig('noisy_line_analysis.png')
plt.close()

print("\nAnalysis of how noise affects the fits:")
print("1. When X has no noise (Case 1), the standard Y on X regression is most appropriate")
print("2. When both X and Y have equal noise (Case 2), both fits are equally valid")
print("3. When X has more noise (Case 3), the reverse fit might be more appropriate")
print("\nRegarding reversibility:")
print("- Theoretically, A' should equal 1/A if the fits were perfectly reversible")
print("- The difference between A' and 1/A increases with noise level")
print("- Higher noise levels generally decrease R² and increase uncertainty in the fits") 