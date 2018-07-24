import matplotlib.pyplot as plt

elements = 6
temperature = [[] for i in range(elements+1)]
cooling_rate = [[] for i in range(elements+1)]
u = []
fu = []
length = 0

for elem in range(elements):
	file_in = open('cooling_output_'+str(elem)+'.dat','r')
	for line in file_in:
		data = line.split()
		temperature[elem+1].append(float(data[0]))
		cooling_rate[elem+1].append(-float(data[1]))
	file_in.close()

#file_in = open('newton_output.dat', 'r')
#for line in file_in:
#	data = line.split()
#	u.append(float(data[0]))
#	fu.append(float(data[1]))
#
#file_in.close()

#p0, = plt.plot(u,fu, color = 'k')
#p1, = plt.plot(u,[0 for x in u], color = 'r')
#p1, = plt.loglog(u,[-x for x in fu], color = 'r')
#plt.xlim([1e13,2.0e14])
#plt.ylim([0e13,2e14])
#plt.xlabel('u')
#plt.ylabel('f(u)')
#plt.show()

#p0, = plt.loglog(temperature[0], cooling_rate[0], linewidth = 0.5, color = 'k', label = 'Total')
p1, = plt.loglog(temperature[1], cooling_rate[1], linewidth = 0.5, color = 'k', label = 'nh = 10^0')
p2, = plt.loglog(temperature[2], cooling_rate[2], linewidth = 0.5, color = 'b', label = 'nh = 10^-1')
p3, = plt.loglog(temperature[3], cooling_rate[3], linewidth = 0.5, color = 'g', label = 'nh = 10^-2')
p4, = plt.loglog(temperature[4], cooling_rate[4], linewidth = 0.5, color = 'r', label = 'nh = 10^-3')
p5, = plt.loglog(temperature[5], cooling_rate[5], linewidth = 0.5, color = 'c', label = 'nh = 10^-4')
p6, = plt.loglog(temperature[6], cooling_rate[6], linewidth = 0.5, color = 'm', label = 'nh = 10^-5')
plt.ylim([1e-24,1e-21])
plt.xlim([1e4,1e8])
plt.legend(handles = [p1,p2,p3,p4,p5,p6])
plt.xlabel('Temperature (K)')
plt.ylabel('Cooling rate (eagle units...)')
plt.show()
