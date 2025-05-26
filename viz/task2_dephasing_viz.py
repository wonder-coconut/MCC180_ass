import numpy as np
import matplotlib.pyplot as plt

def exp_decay(gamma,t):
    return np.exp(-gamma*t)

file = open("../output/dephasing.txt","r").read().split('\n')
params = np.array((file.pop()).split('\t'),dtype='float64')

tau = params[0]
time = params[1]
n_traj = params[2]
gamma_phi = params[3]
delta = params[4]

time_scale = np.arange(0,time,tau)

avg_traj = np.array(file,dtype='float64')

plt.plot(time_scale,avg_traj, label='Numerical average')

#analytical soln
plt.plot(time_scale,exp_decay(2*gamma_phi,time_scale),linestyle='--', color='red',label='Analytical solution')
plt.ylim(-1,1)
plt.xlabel('Time')
plt.ylabel(r'Expectation for $\langle \psi|X|\psi\rangle$')
plt.legend(loc='upper right')
plt.savefig(f'../op_fig/dephase_{int(n_traj)}.png')