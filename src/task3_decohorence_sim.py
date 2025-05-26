import matplotlib.pyplot as plt
import numpy as np
import sys
import lib

#perfectly superposed state
psi = lib.qubit_state_gen(1,1)

gamma_phi = 0.25
gamma = 1
jump_ops = [np.sqrt(gamma_phi) * lib.Z, np.sqrt(gamma) * lib.sigma_minus]

delta = 2
H = delta * lib.Z

tau = 0.01
time = 10

time_scale = np.arange(0,time,tau)

n_traj = int(sys.argv[1])
#trajectories1 = []
trajectories2 = []
#avg_traj1 = 0
avg_traj2 = 0

res = []

#for each wmc process, both output forms (<1|psi> and <psi|X|psi>) were outputted for studying the trends of both under dual jump ops
#the purposes of the assignment only require the latter form, so the former has been commented out
for i in range(n_traj):
    print(i)
    res = lib.wmc_driver(time,psi,tau,jump_ops,H,3).real
    #trajectories1.append(res[:,0])
    trajectories2.append(res[:,1])
    #plt.plot(res[:,0])
    plt.plot(time_scale,res[:,1])

for i in range(n_traj):
    if(i == 0):
        #avg_traj1 = trajectories1[0]
        avg_traj2 = trajectories2[0]
    else:
        #avg_traj1 = avg_traj1 + trajectories1[i]
        avg_traj2 = avg_traj2 + trajectories2[i]

#avg_traj1 /= n_traj
avg_traj2 /= n_traj

#file1 = open("../output/decoherence1.txt","w")
file2 = open("../output/decoherence2.txt","w")

for i in range(len(avg_traj2)):
    #file1.write(str(avg_traj1[i]) + '\n')
    file2.write(str(avg_traj2[i]) + '\n')

file2.write(f'{tau}\t{time}\t{n_traj}\t{gamma}\t{gamma_phi}\t{delta}')

plt.xlabel('Time')
plt.ylabel(r'Expectation for $\langle \psi|X|\psi\rangle$')
plt.savefig(f'../op_fig/decoherence_{n_traj}_unavg.png')
file2.close()
#file1.close()