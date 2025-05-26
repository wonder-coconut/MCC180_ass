import numpy as np
import matplotlib.pyplot as plt

def exp_decay(gamma,t):
    return np.exp(-gamma * t)

#file1 = open("../output/decoherence1.txt","r").read().split('\n')
file2 = open("../output/decoherence2.txt","r").read().split('\n')

#file1.pop()
params = np.array((file2.pop()).split('\t'),dtype='float64')

tau = params[0]
time = params[1]
n_traj = params[2]
gamma = params[3]
gamma_phi = params[4]
delta = params[5]

gamma_res = (2*gamma_phi + gamma/2)
print(gamma_res)

#avg_traj1 = np.array(file1,dtype='float64')
avg_traj2 = np.array(file2,dtype='float64')


time_scale = np.arange(0,10,0.01)
#plt.plot(avg_traj1)
plt.plot(time_scale,avg_traj2, label='Numerical average')

#analytical solution
plt.plot(time_scale,exp_decay(gamma_res,time_scale), label='Analytical solution')
plt.ylim(-1,1)
plt.xlabel('Time')
plt.ylabel(r'Expectation for $\langle \psi|X|\psi\rangle$')
plt.legend(loc='upper right')
plt.savefig(f'../op_fig/decoherence_{int(n_traj)}.png')