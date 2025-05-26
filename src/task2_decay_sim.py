import matplotlib.pyplot as plt
import numpy as np

import lib
import sys

#excited state qubit |1>
psi = lib.qubit_state_gen(0,1)

gamma = 1
time = 10/gamma
tau = 0.01

time_scale = np.arange(0,time,tau)

jump_ops = [np.sqrt(gamma) * lib.sigma_minus]

#number of trajectories
n_traj = int(sys.argv[1])
trajectories = []
avg_traj = 0


for i in range(n_traj):
    print(i)
    res = lib.wmc_driver(time,psi,tau,jump_ops, np.zeros((2,2)),1)
    plt.plot(time_scale,res)
    trajectories.append(res)

for i in range(n_traj):
    if(i == 0):
        avg_traj = trajectories[0]
    else:
        avg_traj = avg_traj + trajectories[i]

avg_traj /= n_traj

file = open("../output/decay.txt","w")

plt.xlabel('Time')
plt.ylabel(r'Population for $\langle 1|\psi\rangle$')
plt.show()

for val in avg_traj:
    file.write(str(val) + '\n')
file.write(f'{tau}\t{time}\t{n_traj}\t{gamma}')
file.close()