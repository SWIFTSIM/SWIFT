import imageio

filename_base = "FastRotor_LR_"

images = []

for ii in range(61):

    filename = filename_base + str(ii).zfill(4) + ".png"
    images.append(imageio.imread(filename))

imageio.mimsave("FastRotor_LR.gif", images)
