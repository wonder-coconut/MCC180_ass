import numpy as np
import matplotlib.pyplot as plt
import cmath

#operator definitions
X = np.array([[0,1],[1,0]],dtype=complex)
Z = np.array([[1,0],[0,-1]],dtype=complex)

H_ad = 1/np.sqrt(2) * np.array([[1,1],[1,-1]],dtype=complex)

#following the convention |psi> = [a,b]T  = a|0> + b|1> 
sigma_minus = np.array([[0,1],[0,0]],dtype=complex)
sigma_plus = np.array([[0,0],[1,0]],dtype=complex)

identity = np.array([[1,0],[0,1]],dtype=complex)

#helper function to return Hermitian conjugate of a vector
def H_conj(v):
    return v.conj().T

#function to generate an arbitrary qubit state: alpha,beta - unnormalized coeffs of |0> and |1> respectively
def qubit_state_gen(alpha=0, beta=0):

    #if no arguments are provided, generate a perfectly superposed state
    if(alpha == 0 and beta == 0):
        alpha = 1/np.sqrt(2)
        beta = 1/np.sqrt(2)
    
    norm = np.sqrt(np.abs(alpha)**2 + np.abs(beta)**2)
    alpha_norm = alpha/norm
    beta_norm = beta/norm

    psi = np.array([[alpha_norm, beta_norm]]).T
    #|psi> = (alpha beta)T
    #taking the convention of ground state being in the 0 position and excited state in 1

    return psi

#evolution operator for mu_bar = 0
def evol_op1(Hamiltonian, jump_ops, tau, p0):
    t1 = identity
    t2 = 1j*tau*Hamiltonian
    t3 = 0
    i = 0
    while(i < len(jump_ops)):
        t3 += np.matmul(H_conj(jump_ops[i]),jump_ops[i])
        i += 1
    t3 = tau/2 * t3
    return(t1 - t2 - t3)/np.sqrt(p0)

#evolution operator for mu_bar != 0
def evol_op2(jump_ops, mu_bar, p):
    numerator = jump_ops[mu_bar-1]
    denominator = np.sqrt(p[mu_bar])
    return numerator/denominator

#function to run a single iteration of the wmc algorithm
def wmc_single_step(psi,tau,jump_ops,H):
    
    #compute probabilites (step 1)
    n = len(jump_ops)

    p = np.zeros(n+1)
    mu = 1
    while(mu <= n):
        left = H_conj(psi) @ H_conj(jump_ops[mu-1])
        right = jump_ops[mu-1] @ psi
        expectation = (left @ right)[0][0].real
        p[mu] = tau * expectation
        mu += 1
    
    p[0] = 1 - np.sum(p)
    
    #compute cumulative probabilites (step 2)
    P = np.zeros(n+2)
    P[0] = 0
    i = 1
    while(i < len(P)):
        P[i] += P[i-1] + p[i-1]
        i+=1
    
    #random nummber (step 3)
    r = np.random.rand()

    #assign mu_bar (step 4)
    mu_bar = 1
    while(mu_bar < len(P)):
        if(P[mu_bar - 1] < r and r < P[mu_bar]):
            break
        mu_bar += 1
    
    #evolve psi (step 5)
    evol_op = 0
    if(mu_bar - 1 == 0):
        evol_op = evol_op1(H,jump_ops,tau,p[0])
    else:
        evol_op = evol_op2(jump_ops, mu_bar-1, p)
    
    psi_evol = evol_op @ psi
    psi_evol /= np.linalg.norm(psi_evol)

    return psi_evol

#function to run all the time steps of the wmc algorithm
def wmc_driver(time, psi, tau, jump_ops, H, op_flag):
    
    iterations = int(time/tau)
    res = []
    for i in range(iterations):
        psi_evol = wmc_single_step(psi, tau, jump_ops, H)
        if(op_flag == 1):
            #decay output
            res.append(np.abs(psi_evol[1, 0])**2) #measuring (<1|psi>)^2
        elif(op_flag == 2):
            #dephasing output
            res.append((H_conj(psi_evol) @ X @ psi_evol)[0][0].real) #measuring <psi|X|psi>
        elif(op_flag == 3):
            #decoherence output
            temp = np.array([np.abs(psi_evol[1, 0])**2, (H_conj(psi_evol) @ X @ psi_evol)[0][0]])
            res.append(temp)
        psi = psi_evol
    
    
    return np.array(res)