import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

n_snaps = 100
directory = '../'
snap_name_base = "uniformBox_"
plot_directory = "./"
plot_filename_base = "snap"
for i in range(n_snaps):
    snap_number = "%03d" %i
    filename = directory + snap_name_base + snap_number + ".hdf5"
    f = h5.File(filename)
    #find the box size

    header = f["Header"]
    box_size = header.attrs["BoxSize"][0]
    pos = np.array(f["PartType0/Coordinates"])
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    #x_projection
    plt.plot(y,z,'bo')
    plt.xlabel("y")
    plt.ylabel("z")
    plt.xlim((0.0,box_size))
    plt.ylim((0.0,box_size))
    plot_dir = "%s/projection_x" %plot_directory
    plot_filename = "%s/%s_%03d.png" %(plot_dir,plot_filename_base,i)
    plt.savefig(plot_filename,format = "png")
    plt.close()

    #y_projection
    plt.plot(x,z,'bo')
    plt.xlabel("x")
    plt.ylabel("z")
    plt.xlim((0.0,box_size))
    plt.ylim((0.0,box_size))
    plot_dir = "%s/projection_y" %plot_directory
    plot_filename = "%s/%s_%03d.png" %(plot_dir,plot_filename_base,i)
    plt.savefig(plot_filename,format = "png")
    plt.close()

    #x_projection
    plt.plot(x,y,'bo')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xlim((0.0,box_size))
    plt.ylim((0.0,box_size))
    plot_dir = "%s/projection_z" %plot_directory
    plot_filename = "%s/%s_%03d.png" %(plot_dir,plot_filename_base,i)
    plt.savefig(plot_filename,format = "png")
    plt.close()
    

