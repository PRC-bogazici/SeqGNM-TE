import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Enter Inputs
name_EC = 'b03.xlsx'  # Name of Excel file containing evolutionary contacts
N = 360  # Number of residues
m1 = 1  # First mode
m2 = 10  # Last mode
tau_0 = 6.0
threshold = 0.80

# Prepare Kirchhoff matrix
# Load Excel file
df_EC = pd.read_excel(name_EC, header=None)

# Coerce all columns to numeric floats
df_EC = df_EC.apply(pd.to_numeric, errors='coerce')
EC = df_EC.dropna().to_numpy()

K = np.zeros((N, N))  # Kirchhoff matrix

# Coupled contacts
for row in EC:
    if row[3] > threshold:
        i_idx = int(row[0]) - 1
        j_idx = int(row[1]) - 1
        val = float(row[2])
        K[i_idx, j_idx] = -val
        K[j_idx, i_idx] = -val

# Chain connectivity
for i in range(N):
    for j in range(N):
        if i != j and abs(i - j) < 4:
            K[i, j] -= 1

# Diagonal elements (-sum of non-diagonal elements)
np.fill_diagonal(K, -np.sum(K, axis=1))

# Matrix Decomposition (SVD)
U, S, Vt = np.linalg.svd(K)

# Sort eigenvalues and eigenvectors in ascending order
S = S[::-1]
U = U[:, ::-1]

# Find Maximum Information Tau
lower = 0.01
inc = 20.0
upper = 4000.0
tau_range = np.arange(lower, upper + inc, inc)

Tau = np.zeros((N, N))

i_residues = np.arange(1, N, 10)
j_residues = np.arange(4, N, 10)

# Mode bounds
k_start = m1
k_end = m2 + 1  #
k_indices = np.arange(k_start, k_end)

for i in i_residues:
    for j in j_residues:
        Tij = np.zeros(len(tau_range))

        for idx, tau in enumerate(tau_range):
            u_j = U[j, k_indices]
            u_i = U[i, k_indices]
            s_k = S[k_indices]

            exp_term = np.exp(-s_k * tau / tau_0)

            sum1 = np.sum(u_j * (1.0 / s_k) * u_j)
            sum2 = np.sum(u_j * (1.0 / s_k) * u_j * exp_term)
            sum3 = np.sum(u_i * (1.0 / s_k) * u_i)
            sum4 = np.sum(u_i * (1.0 / s_k) * u_j)
            sum5 = np.sum(u_i * (1.0 / s_k) * u_j * exp_term)

            term1 = sum1 ** 2 - sum2 ** 2
            term2 = sum3 * sum1 ** 2 + 2 * sum4 * sum2 * sum5 - (sum5 ** 2 + sum4 ** 2) * sum1 - sum2 ** 2 * sum3
            term3 = sum1
            term4 = sum3 * sum1 - sum4 ** 2

            # Avoid computational crash / negative values inside log
            if abs(term1) < 9e-4 and term1 < 0:
                term1 = abs(term1)
            if abs(term2) < 9e-4 and term2 < 0:
                term2 = abs(term2)
            if abs(term3) < 9e-4 and term3 < 0:
                term3 = abs(term3)
            if abs(term4) < 9e-4 and term4 < 0:
                term4 = abs(term4)

            term1 = 0.5 * np.log(term1)
            term2 = 0.5 * np.log(term2)
            term3 = 0.5 * np.log(term3)
            term4 = 0.5 * np.log(term4)

            Tij[idx] = term1 - term2 - term3 + term4

        # Detect tau giving max information transfer
        peaks, _ = find_peaks(Tij)
        if len(peaks) > 0:
            Tau[i, j] = tau_range[peaks[0]]

taus = Tau[Tau != 0]
tau = np.mean(taus) * 3.0  # Optimal time delay

# Calculate Transfer Entropy
Tij = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        if i == j:
            Tij[i, j] = 0.0
        else:
            u_j = U[j, k_indices]
            u_i = U[i, k_indices]
            s_k = S[k_indices]

            exp_term = np.exp(-s_k * tau / tau_0)

            sum1 = np.sum(u_j * (1.0 / s_k) * u_j)
            sum2 = np.sum(u_j * (1.0 / s_k) * u_j * exp_term)
            sum3 = np.sum(u_i * (1.0 / s_k) * u_i)
            sum4 = np.sum(u_i * (1.0 / s_k) * u_j)
            sum5 = np.sum(u_i * (1.0 / s_k) * u_j * exp_term)

            term1 = sum1 ** 2 - sum2 ** 2
            term2 = sum3 * sum1 ** 2 + 2 * sum4 * sum2 * sum5 - (sum5 ** 2 + sum4 ** 2) * sum1 - sum2 ** 2 * sum3
            term3 = sum1
            term4 = sum3 * sum1 - sum4 ** 2

            if abs(term1) < 9e-4 and term1 < 0:
                term1 = abs(term1)
            if abs(term2) < 9e-4 and term2 < 0:
                term2 = abs(term2)
            if abs(term3) < 9e-4 and term3 < 0:
                term3 = abs(term3)
            if abs(term4) < 9e-4 and term4 < 0:
                term4 = abs(term4)

            term1 = 0.5 * np.log(term1)
            term2 = 0.5 * np.log(term2)
            term3 = 0.5 * np.log(term3)
            term4 = 0.5 * np.log(term4)

            Tij[i, j] = term1 - term2 - term3 + term4

# %% Calculate Net Transfer Entropy
netTE = Tij - Tij.T

# %% Figure Plotting
plt.figure(figsize=(7, 6))
im = plt.imshow(netTE, origin='lower', cmap='jet', vmin=-0.01, vmax=0.01)

plt.colorbar(im)
plt.xlabel('Residue Index (Affected)', fontsize=14)
plt.ylabel('Residue Index (Effector)', fontsize=14)

if m1 != m2:
    plot_name = f"{m1}-{m2} Modes - Net Transfer Entropy"
else:
    plot_name = f"Mode {m1} - Net Transfer Entropy"

plt.title(plot_name, fontsize=14)
plt.tight_layout()
plt.show()

# Find Collectivity
net_pos = np.maximum(netTE, 0)
net_sq = net_pos ** 2
sum_sq = np.sum(net_sq, axis=1, keepdims=True)

# Avoid division by zero for rows with zero sum
with np.errstate(divide='ignore', invalid='ignore'):
    alfa1 = np.where(sum_sq > 0, 1.0 / sum_sq, 0.0)

P = net_sq * alfa1

P_logP = np.zeros_like(P)
pos_mask = P > 0
P_logP[pos_mask] = P[pos_mask] * np.log(P[pos_mask])

insidesum = np.sum(P_logP, axis=1)
col = np.exp(-insidesum) / N  # Collectivity vector