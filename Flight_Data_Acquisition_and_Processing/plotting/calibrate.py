import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------
# INPUT DATA (EDIT THESE)
# ----------------------------------------------------
# Measured values (sensor output)
x_meas = np.array([
    101325.0,
    101300.2,
    101275.6,
    101250.1,
    101225.0,
    101200.3
])

# True / nominal values (ground truth)
y_true = np.array([
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5
])

# Polynomial degree
DEGREE = 5

# ----------------------------------------------------
# FIT POLYNOMIAL
# ----------------------------------------------------
coeffs = np.polyfit(x_meas, y_true, DEGREE)

# Polynomial function
poly = np.poly1d(coeffs)

# Predicted values
y_fit = poly(x_meas)

# ----------------------------------------------------
# ERROR METRICS
# ----------------------------------------------------
error = y_fit - y_true
rmse = np.sqrt(np.mean(error**2))
max_err = np.max(np.abs(error))

# ----------------------------------------------------
# PRINT RESULTS
# ----------------------------------------------------
print("\n=== Polynomial Calibration (degree {}) ===".format(DEGREE))
print("y = a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0\n")

for i, c in enumerate(coeffs):
    power = DEGREE - i
    print(f"a{power} = {c:.10e}")

print("\n--- Error metrics ---")
print(f"RMSE     = {rmse:.6f}")
print(f"Max err  = {max_err:.6f}")

# ----------------------------------------------------
# OPTIONAL: PLOT
# ----------------------------------------------------
x_dense = np.linspace(min(x_meas), max(x_meas), 500)
y_dense = poly(x_dense)

plt.figure(figsize=(7,4))
plt.plot(x_meas, y_true, 'o', label="Ground truth")
plt.plot(x_dense, y_dense, '-', label="Polynomial fit")
plt.xlabel("Measured value")
plt.ylabel("True value")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
