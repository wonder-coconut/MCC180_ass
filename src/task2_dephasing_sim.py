import matplotlib.pyplot as plt
import numpy as np
import sys
import lib

#perfectly superposed state
psi = lib.qubit_state_gen(1,1)

gamma_phi = 0.25
jump_ops = [np.sqrt(gamma_phi) * lib.Z]

#hamiltonian
delta = 1.5
H = delta * lib.Z

tau = 0.01
time = 10
time_scale = np.arange(0,time,tau)

n_traj = int(sys.argv[1])
trajectories = []
avg_traj = 0

for i in range(n_traj):
    print(i)
    res = lib.wmc_driver(time, psi, tau, jump_ops, H,2)
    plt.plot(time_scale,res)
    trajectories.append(res)

for i in range(n_traj):
    if(i == 0):
        avg_traj = trajectories[0]
    else:
        avg_traj = avg_traj + trajectories[i]

avg_traj /= n_traj

plt.xlabel('Time')
plt.ylabel(r'Expectation for $\langle \psi|X|\psi\rangle$')
plt.savefig(f'../op_fig/dephase_{n_traj}_unavg.png')

file = open("../output/dephasing.txt","w")

for val in avg_traj:
    file.write(str(val) + '\n')
file.write(f'{tau}\t{time}\t{n_traj}\t{gamma_phi}\t{delta}')
file.close()