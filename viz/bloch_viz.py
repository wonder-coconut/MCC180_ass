import qutip 
import numpy as np
import matplotlib.pyplot as plt

b = qutip.Bloch()

psi_0 = 1/np.sqrt(2) * (qutip.basis(2, 0) + qutip.basis(2,1))
psi_1 = 1/np.sqrt(2) * (qutip.basis(2, 0) - qutip.basis(2,1)) 
b.add_states([psi_0,psi_1])

b.render()
plt.savefig('../op_fig/bloch.png')