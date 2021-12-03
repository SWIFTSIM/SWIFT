from swiftsimio import load
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.animation import FuncAnimation

import os

img_files = [for f in os.listdir(os.path.abspath(__file__)) if (f[-4:] == 'hdf5' and f != 'uniformBox.hdf5')]

print(img_files)

quit()

# redraw_fn draw frame f
def redraw_fn(f, axes):
    img_file = img_files[f]
    
    data = load(img_file)

    x_gas = data.gas.coordinates
    v_gas = data.gas.velocities
    P = data.gas.pressures

	 rho_gas = data.gas.densities	
	 B_gas = data.gas.magnetic_flux_density
	 
    if not redraw_fn.initialized:    	
      ax = plt.figure().add_subplot(projection='3d') 
		ax.quiver(x_gas[:, 0],x_gas[:, 1],x_gas[:, 2],B_gas[:, 0],B_gas[:, 1],B_gas[:, 2], length=0.05, normalize=True)
		ax.set_xlabel('x')
      ax.set_ylabel('y')
      plt.show()
      redraw_fn.initialized = True
    else:
    	ax = plt.figure().add_subplot(projection='3d') 
		ax.quiver(x_gas[:, 0],x_gas[:, 1],x_gas[:, 2],B_gas[:, 0],B_gas[:, 1],B_gas[:, 2], length=0.05, normalize=True)
		ax.set_xlabel('x')
      ax.set_ylabel('y')
      plt.show()
      redraw_fn.im.set_array(img)
redraw_fn.initialized = False

videofig(len(img_files), redraw_fn, play_fps=30)

quit()
