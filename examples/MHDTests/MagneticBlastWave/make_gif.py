import imageio

filename_base = "MagneticBlastWave_LR_"

images = []

for ii in range(41):

	filename = filename_base + str(ii).zfill(4) + ".png"
	images.append(imageio.imread(filename))

imageio.mimsave('MagneticBlastWave_LR.gif', images)
