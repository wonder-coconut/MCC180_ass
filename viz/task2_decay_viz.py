import numpy as np
import matplotlib.pyplot as plt

def exp_decay(gamma,t):
    return np.exp(-gamma*t)

file = open("../output/decay.txt","r").read().split('\n')
params = np.array((file.pop()).split('\t'),dtype='float64')

tau = params[0]
time = params[1]
n_traj = params[2]
gamma = params[3]

time_scale = np.arange(0,time,tau)

avg_traj = np.array(file,dtype='float64')

#finding jump timestamp
#for i in range(len(avg_traj)):
#    if(avg_traj[i+1] - avg_traj[i] == -1):
#        break
#print(i)

plt.plot(time_scale,avg_traj,label='Averaged Trajectory')

#analytical soln
plt.plot(time_scale, exp_decay(gamma,time_scale), linestyle='--', color='red', label='Analytical Solution')


plt.xlabel('Time')
plt.ylabel(r'Population for $\langle 1|\psi\rangle$')
plt.legend(loc='upper right')
plt.savefig(f'../op_fig/decay_{int(n_traj)}.png')