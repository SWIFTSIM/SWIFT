from swiftsimio import load
import numpy as np
import matplotlib.pyplot as plt

# Of course, replace this path with your own snapshot should you be using
# custom data.
data = load("uniformBox_0008.hdf5")

# Note that gas_properties is an instance of SWIFTParticleTypeMetadata
print(data.metadata.gas_properties.field_names)

x_gas = data.gas.coordinates
v_gas = data.gas.velocities
P = data.gas.pressures

rho_gas = data.gas.densities	
rhosq_gas = data.gas.densities_squared
B_gas = data.gas.magnetic_flux_density

"""
plt.figure()
plt.subplot(211)
plt.scatter(x_gas[:, 0], np.sqrt(rhosq_gas)/rho_gas, s = 1., c='r')
plt.xlim([.5,1.5])
plt.ylim([0.95,1.05])
plt.axhline(y=1., color='k', linestyle='-')
plt.ylabel(r'$\sqrt{\widetilde{\rho^2}} / \widetilde{\rho}$')

plt.subplot(212)
plt.scatter(x_gas[:,0], rho_gas, s = 1., c='r')
plt.xlim([.5,1.5])
plt.xlabel(r'$x$')
plt.ylabel(r'$\bar{\rho}$')
plt.show()
"""

ax = plt.figure().add_subplot(projection='3d') 
ax.quiver(x_gas[:, 0],x_gas[:, 1],x_gas[:, 2],B_gas[:, 0],B_gas[:, 1],B_gas[:, 2], length=0.05, normalize=True)
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
