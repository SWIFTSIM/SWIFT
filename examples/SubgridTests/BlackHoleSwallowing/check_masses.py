import h5py as h5
import numpy as np


f = h5.File("bh_swallowing_serial_0001.hdf5", "r")
f_mpi = h5.File("bh_swallowing_mpi_0001.hdf5", "r")
f_ics = h5.File("bh_swallowing.hdf5", "r")

m = np.array(f["/PartType5/Masses"])#[:]
m_mpi = np.array(f_mpi["/PartType5/Masses"])#[:]
m_ics = np.array(f_ics["/PartType5/Masses"])#[:]

ids = np.array(f["/PartType5/ParticleIDs"])#[:]
ids_mpi = np.array(f_mpi["/PartType5/ParticleIDs"])#[:]
ids_ics = np.array(f_ics["/PartType5/ParticleIDs"])#[:]

rank = np.argsort(ids)
rank_mpi = np.argsort(ids_mpi)
rank_ics = np.argsort(ids_ics)

m = m[rank]
m_mpi = m_mpi[rank_mpi]
m_ics = m_ics[rank_ics]

ids = ids[rank]
ids_mpi = ids_mpi[rank_mpi]
ids_ics = ids_ics[rank_ics]


# Check
#print ids - ids_mpi
#print ids - ids_ics

count_wrong = 0

print "ID -- ICs mass -- serial mass -- MPI mass"
for i in range(np.size(ids)):
    print "%14d %5f %5f %5f"%(ids[i], m_ics[i], m[i], m_mpi[i]),
    
    if m[i] > m_ics[i] * 1.1:
        print "#",
    else:
        print " ",
    
    if m[i] != m_mpi[i]:
        count_wrong += 1
        print " <---"
    else:
        print ""

print "Found", count_wrong, "wrong BH masses."


ids_gas_mpi = f_mpi["/PartType0/ParticleIDs"][:]
ids_removed_mpi = np.loadtxt("list_mpi.txt")

for i in range(np.size(ids_removed_mpi)):
    result = np.where(ids_gas_mpi == ids_removed_mpi)
    print result


ids_gas = f["/PartType0/ParticleIDs"][:]
ids_removed = np.loadtxt("list.txt")

for i in range(np.size(ids_removed)):
    result = np.where(ids_gas == ids_removed)
    print result

#rho_gas = f["/PartType0/Densities"][:]
#print np.mean(rho_gas), np.std(rho_gas)
