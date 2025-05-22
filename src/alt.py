import numpy as np
import matplotlib.pyplot as plt

# --- Constants and operators ---
# Basis states
zero = np.array([[1], [0]], dtype=complex)
one = np.array([[0], [1]], dtype=complex)

# Pauli matrices
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Hadamard operator
H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

# Annihilation (sigma-) and creation (sigma+) operators
sigma_minus = np.array([[0, 1], [0, 0]], dtype=complex)
sigma_plus = np.array([[0, 0], [1, 0]], dtype=complex)

# --- State initialization ---
def create_qubit_state(alpha, beta):
    state = np.array([[alpha], [beta]], dtype=complex)
    norm = np.linalg.norm(state)
    return state / norm

# --- Inner product helper ---
def inner_product(state1, op, state2):
    return (state1.conj().T @ op @ state2).item()

# --- WFMC: single time step ---
def wfmc_step(psi, H, collapse_ops, tau):
    probs = [tau * inner_product(psi, c.conj().T @ c, psi).real for c in collapse_ops]
    p0 = 1 - sum(probs)
    P = [p0]
    for p in probs:
        P.append(P[-1] + p)
    
    r = np.random.rand()
    mu_bar = next(i for i, val in enumerate(P) if r < val) - 1
    
    if mu_bar == -1:  # No jump
        effective_H = H - 0.5j * sum(c.conj().T @ c for c in collapse_ops)
        psi_new = (np.eye(2) - 1j * tau * effective_H) @ psi
        psi_new /= np.linalg.norm(psi_new)
    else:  # Jump occurred
        c_mu = collapse_ops[mu_bar]
        psi_new = c_mu @ psi
        psi_new /= np.linalg.norm(psi_new)
    
    return psi_new

# --- WFMC: full trajectory ---
def simulate_trajectory(psi0, H, collapse_ops, tau, t_max):
    steps = int(t_max / tau)
    psi = psi0.copy()
    probs_1 = [np.abs(one.conj().T @ psi)**2]

    for _ in range(steps):
        psi = wfmc_step(psi, H, collapse_ops, tau)
        probs_1.append(np.abs(one.conj().T @ psi)**2)

    return probs_1


# Parameters
gamma = 1.0
tau = 0.01
t_max = 10 / gamma
H = np.zeros((2, 2), dtype=complex)  # No Hamiltonian for decay-only case
collapse_ops = [np.sqrt(gamma) * sigma_minus]

# Initial state: |1>
psi0 = zero.copy()

t_arr = []
for i in range(100):
    trajectory = simulate_trajectory(psi0, H, collapse_ops, tau, t_max)
    trajectory = np.array(trajectory)[:,0,0]
    t_arr.append(trajectory)

for i in range(100):
    plt.plot(t_arr[i])

plt.show()