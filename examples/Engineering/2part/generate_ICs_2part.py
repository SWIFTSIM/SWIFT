from swiftsimio import Writer
import swiftsimio.metadata.particle as swp
import swiftsimio.metadata.writer.required_fields as swmw
import swiftsimio.metadata.unit.unit_fields as swuf
import numpy as np
import unyt
import pdb


def generate_units(m, l, t, I, T):
    dict_out = swuf.generate_units(m , l, t, I, T)

    boundary = {
            "coordinates" : l,
            "masses": m,
            "particle_ids": None,
            "velocities": l / t, 
            "potential": l * l / (t * t),
            "density": m / (l ** 3),
            "entropy": m * l ** 2 / (t ** 2 * T),
            "internal_energy": (l / t) ** 2,
            "smoothing_length": l,
            "pressure": m / (l * t ** 2),
            "diffusion": None,
            "sfr": m / t,
            "temperature": T,
            "viscosity": None,
            "specific_sfr": 1 / t,
            "material_id": None,
            "diffusion": None,
            "viscosity": None,
            "radiated_energy": m * (l / t) ** 2,
            "constant_acceleration": l / (t*t),
            "is_boundary": None,
            }
    fluid = { **boundary }

    dict_out["boundary"] = boundary
    dict_out["fluid"] = fluid
    return dict_out


#Use default units, i.e. cm, grams, seconds, Ampere, Kelvin
unit_system = unyt.UnitSystem(name="default", length_unit=unyt.m, mass_unit=unyt.kg, time_unit=unyt.s)

#Add boundary particle type
swp.particle_name_underscores[6] = "boundary"
swp.particle_name_class[6] = "Boundary"
swp.particle_name_text[6] = "Boundary"

swp.particle_name_underscores[7] = "fluid"
swp.particle_name_class[7] = "Fluid"
swp.particle_name_text[7] = "Fluid"

swp.particle_fields.boundary = { 'constant_acceleration' : 'ConstantAcceleration',
                                  'is_boundary' : 'IsBoundary',
                                   **swp.particle_fields.gas
                               }
swp.particle_fields.fluid = { **swp.particle_fields.boundary }


swmw.fluid = {'smoothing_length' : 'SmoothingLength', 
              'constant_acceleration' : 'ConstantAcceleration',
              'is_boundary' : 'IsBoundary', 'density': 'Density' ,**swmw.gas, **swmw.shared}
swmw.boundary = { **swmw.fluid }



n_hydro = 2
n_p = n_hydro 

boxsize = 1.5*unyt.m

x = Writer(unit_system, boxsize, unit_fields_generate_units=generate_units)



x.fluid.coordinates = np.zeros((n_hydro,3)) * unyt.m
x.fluid.coordinates[0][0] = 0.5 * unyt.m
x.fluid.coordinates[0][1] = 0.25 * unyt.m
x.fluid.coordinates[0][2] = 0.5 * unyt.m
x.fluid.coordinates[1][0] = 0.5 * unyt.m
x.fluid.coordinates[1][1] = 1.25 * unyt.m
x.fluid.coordinates[1][2] = 0.5 * unyt.m


x.fluid.velocities = np.zeros((n_hydro,3)) * unyt.m/unyt.s

x.fluid.masses = np.ones(n_hydro, dtype = float) * 10.0* unyt.kg
x.fluid.density = np.ones(n_hydro, dtype = float) * (1000 * unyt.kg / (unyt.m*unyt.m*unyt.m))
x.fluid.velocities[0][1] = 1.0 * unyt.m / unyt.s
x.fluid.velocities[1][1] = -1.0 * unyt.m / unyt.s

x.fluid.smoothing_length = np.ones(n_hydro, dtype = float) * (0.13 * unyt.m)

x.fluid.is_boundary = np.zeros(n_hydro, dtype = int)
x.fluid.constant_acceleration = np.zeros((n_hydro,3)) * unyt.m/(unyt.s*unyt.s)
x.fluid.internal_energy = np.zeros((n_hydro)) * (unyt.m/unyt.s)*(unyt.m/unyt.s)


#x.gas.coordinates = np.zeros((n_hydro,3)) * unyt.cm
#x.gas.coordinates[0][0] = 50.0 * unyt.cm
#x.gas.coordinates[0][1] = 75.0 * unyt.cm
#x.gas.coordinates[0][2] = 50.0 * unyt.cm
#
#
#x.gas.velocities = np.zeros((n_hydro,3)) * unyt.cm/unyt.s
#
#x.gas.masses = np.ones(n_hydro, dtype = float) * unyt.g
#
#x.gas.smoothing_length = np.ones(n_hydro, dtype = float) * (5.0 * unyt.cm)
#x.gas.internal_energy =  np.ones(n_hydro, dtype=float) * (293 * unyt.J/unyt.g)

x.write("2part_test.hdf5")
