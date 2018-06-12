import matplotlib.pyplot as plt

elements = 11
temperature = [[] for i in range(10)]
cooling_rate = [[] for i in range(10)]
heating_rate = [[] for i in range(10)]
dlambda_du = [[] for i in range(10)]
u = []
fu = []
length = 0
for i in range(10):
        file_in = open('cooling_integration_output_'+str(i)+'.dat','r')
	for line in file_in:
		data = line.split()
		temperature[i].append(float(data[0]))
		cooling_rate[i].append(-float(data[1]))
		heating_rate[i].append(float(data[1]))
		dlambda_du[i].append(-float(data[2]))

file_in.close()

#p10, = plt.loglog(temperature[0], dlambda_du[0], color = 'k', 		label = 'Heating')
#p11, = plt.loglog(temperature[1], dlambda_du[1], color = 'r', 		label = 'Heating')
#p12, = plt.loglog(temperature[2], dlambda_du[2], color = 'g', 		label = 'Heating')
#p13, = plt.loglog(temperature[3], dlambda_du[3], color = 'b', 		label = 'Heating')
#p14, = plt.loglog(temperature[4], dlambda_du[4], color = 'c', 		label = 'Heating')
#p15, = plt.loglog(temperature[5], dlambda_du[5], color = 'm', 		label = 'Heating')
#p16, = plt.loglog(temperature[6], dlambda_du[6], color = 'y', 		label = 'Heating')
#p17, = plt.loglog(temperature[7], dlambda_du[7], color = 'olive', 	label = 'Heating')
#p18, = plt.loglog(temperature[8], dlambda_du[8], color = 'orchid', 	label = 'Heating')
#p19, = plt.loglog(temperature[9], dlambda_du[9], color = 'orange', 	label = 'Heating')
##plt.ylim([-1e-24,1e-24])
#plt.xlim([1e2,1e9])
#plt.xlabel('Temperature (K)')
#plt.ylabel('dLambda_du')
#plt.show()
##plt.savefig("dlambda_du.png")

p0,  = plt.loglog(temperature[0], cooling_rate[0], 		     color = 'k', 		label = 'Cooling')
#p1,  = plt.loglog(temperature[1], cooling_rate[1], 		     color = 'r', 		label = 'Cooling')
#p2,  = plt.loglog(temperature[2], cooling_rate[2], 		     color = 'g', 		label = 'Cooling')
#p3,  = plt.loglog(temperature[3], cooling_rate[3], 		     color = 'b', 		label = 'Cooling')
#p4,  = plt.loglog(temperature[4], cooling_rate[4], 		     color = 'c', 		label = 'Cooling')
#p5,  = plt.loglog(temperature[5], cooling_rate[5], 		     color = 'm', 		label = 'Cooling')
#p6,  = plt.loglog(temperature[6], cooling_rate[6], 		     color = 'y', 		label = 'Cooling')
#p7,  = plt.loglog(temperature[7], cooling_rate[7], 		     color = 'olive', 		label = 'Cooling')
#p8,  = plt.loglog(temperature[8], cooling_rate[8], 		     color = 'orchid',		label = 'Cooling')
#p9,  = plt.loglog(temperature[9], cooling_rate[9], 		     color = 'orange',		label = 'Cooling')
p10, = plt.loglog(temperature[0], heating_rate[0], linestyle = '--', color = 'k', 		label = 'Heating')
#p11, = plt.loglog(temperature[1], heating_rate[1], linestyle = '--', color = 'r', 		label = 'Heating')
#p12, = plt.loglog(temperature[2], heating_rate[2], linestyle = '--', color = 'g', 		label = 'Heating')
#p13, = plt.loglog(temperature[3], heating_rate[3], linestyle = '--', color = 'b', 		label = 'Heating')
#p14, = plt.loglog(temperature[4], heating_rate[4], linestyle = '--', color = 'c', 		label = 'Heating')
#p15, = plt.loglog(temperature[5], heating_rate[5], linestyle = '--', color = 'm', 		label = 'Heating')
#p16, = plt.loglog(temperature[6], heating_rate[6], linestyle = '--', color = 'y', 		label = 'Heating')
#p17, = plt.loglog(temperature[7], heating_rate[7], linestyle = '--', color = 'olive', 		label = 'Heating')
#p18, = plt.loglog(temperature[8], heating_rate[8], linestyle = '--', color = 'orchid', 		label = 'Heating')
#p19, = plt.loglog(temperature[9], heating_rate[9], linestyle = '--', color = 'orange', 		label = 'Heating')
plt.ylim([1e7,1e20])
plt.xlim([1e10,1e19])
plt.xlabel('Temperature (K)')
plt.ylabel('f(x)')
##plt.legend(handles = [p0,p1])
plt.show()
#plt.savefig("fu_1.png")
