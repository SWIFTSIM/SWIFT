import matplotlib.pyplot as plt

elements = 9
temperature = []
cooling_rate = [[] for i in range(elements+1)]
length = 0
file_in = open('cooling_output.dat', 'r')
for line in file_in:
	data = line.split()
	temperature.append(float(data[0]))
	cooling_rate[0].append(-float(data[1]))

file_in.close()

for elem in range(elements):
	file_in = open('cooling_output_'+str(elem)+'.dat','r')
	for line in file_in:
        	data = line.split()
        	cooling_rate[elem+1].append(-float(data[0]))
	file_in.close()


p0, = plt.loglog(temperature, cooling_rate[0], linewidth = 0.5, color = 'k', label = 'Total')
p1, = plt.loglog(temperature, cooling_rate[1], linewidth = 0.5, color = 'b', label = 'Carbon')
p2, = plt.loglog(temperature, cooling_rate[2], linewidth = 0.5, color = 'g', label = 'Nitrogen')
p3, = plt.loglog(temperature, cooling_rate[3], linewidth = 0.5, color = 'r', label = 'Oxygen')
p4, = plt.loglog(temperature, cooling_rate[4], linewidth = 0.5, color = 'c', label = 'Neon')
p5, = plt.loglog(temperature, cooling_rate[5], linewidth = 0.5, color = 'm', label = 'Magnesium')
p6, = plt.loglog(temperature, cooling_rate[6], linewidth = 0.5, color = 'y', label = 'Silicon')
p7, = plt.loglog(temperature, cooling_rate[7], linewidth = 0.5, color = 'lightgray', label = 'Sulphur')
p8, = plt.loglog(temperature, cooling_rate[8], linewidth = 0.5, color = 'olive', label = 'Calcium')
p9, = plt.loglog(temperature, cooling_rate[9], linewidth = 0.5, color = 'saddlebrown', label = 'Iron')
plt.ylim([1e-24,1e-21])
plt.xlim([1e4,1e8])
plt.xlabel('Temperature (K)')
plt.ylabel('Cooling rate (eagle units...)')
plt.legend(handles = [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9])
plt.show()
