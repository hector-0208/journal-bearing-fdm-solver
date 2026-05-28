import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson


def analyze_journal_bearing(epsilon, L_D_ratio, M=88, N=22):
    D_L_ratio = 1.0 / L_D_ratio
    theta_0, theta_n = 0.0, 2.0 * np.pi
    z_0, z_n = 0.0, 1.0
    # Calculate step sizes
    d_theta = (theta_n - theta_0) / (M - 1)
    dz = (z_n - z_0) / (N - 1)
    # Create grid arrays
    theta = np.linspace(theta_0, theta_n, M)
    z_bar = np.linspace(z_0, z_n, N)
    # Pre-compute the constant multiplier 'C'
    C = 0.25 * (D_L_ratio ** 2) * ((d_theta / dz) ** 2)
    denominator = 2.0 * (1.0 + C)
    # Initialize pressure array with boundary conditions (zeros at edges)
    P = np.zeros((M, N))

    max_iteration = 3000
    tolerance = 1e-5

    for iteration in range(max_iteration):
        max_error = 0.0
        
        for i in range(1, M - 1):
            # Film thickness equation: h_bar = 1 + epsilon * cos(theta)
            h_i = 1.0 + epsilon * np.cos(theta[i])
            for j in range(1, N - 1):
                P_old = P[i, j]
                term1 = P[i+1, j] + P[i-1, j]
                term2 = - (3.0 * d_theta / (2.0 * h_i)) * epsilon * np.sin(theta[i]) * (P[i+1, j] - P[i-1, j])
                term3 = (epsilon * np.sin(theta[i]) / (h_i ** 3)) * (d_theta ** 2)
                term4 = C * (P[i, j+1] + P[i, j-1])
                P_new = (term1 + term2 + term3 + term4) / denominator
                # Apply Reynolds Cavitation Boundary Condition
                if P_new < 0:
                    P_new = 0.0
                    
                P[i, j] = P_new
                error = abs(P_new - P_old)
                max_error = max(max_error, error)
        if max_error < tolerance:
            break

    # Numerical Integration (Simpson's 1/3rd Rule) ---
    P_cos = P * np.cos(theta)[:, np.newaxis]
    P_sin = P * np.sin(theta)[:, np.newaxis]

    # Integrate first along the axial direction (z) for all theta nodes
    Wy_theta = simpson(P_cos, x=z_bar, axis=1)
    Wx_theta = simpson(P_sin, x=z_bar, axis=1)

    # Integrate the resulting 1D arrays along the circumferential direction (theta)
    W_y = simpson(Wy_theta, x=theta)
    W_x = simpson(Wx_theta, x=theta)
    
    W_y = -W_y # Negative because pressure acts inward

    # Calculate total non-dimensional load capacity and attitude angle
    W_bar = np.sqrt(W_x**2 + W_y**2)
    phi = np.degrees(np.arctan2(W_x, W_y))

    return P, theta, W_bar, phi

# Generate 2D Plot for Specific Case
epsilon_plot = 0.9
LD_plot = 0.5
P_plot, theta_plot, _, _ = analyze_journal_bearing(epsilon_plot, LD_plot)

center_index = P_plot.shape[1] // 2

plt.figure(figsize=(8, 5))
plt.plot(theta_plot, P_plot[:, center_index], linewidth=2, label="Pressure at center (z=0.5)")
plt.title(f"Pressure Distribution at Bearing Center (L/D = {LD_plot}, e = {epsilon_plot})")
plt.xlabel("Circumferential Angle (theta)[rad]")
plt.ylabel("Non-Dimensional Pressure")
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()

plt.show()

# Generate Tables for Multiple Configurations
LD_ratios = [0.5, 1.0, 2.0]
eccentricities = [0.1, 0.3, 0.5, 0.7, 0.9]

for LD in LD_ratios:
    print(f"\nTABLE FOR L/D = {LD}\n")
    print(f"{'Eccentricity (e)':<20} | {'Load Capacity (W_bar)':<20} | {'Attitude Angle (deg)':<20}")
    for e in eccentricities:
        _, _, W, phi = analyze_journal_bearing(e, LD)
        print(f"{e:<20.1f} | {W:<20.4f} | {phi:<20.2f}")
