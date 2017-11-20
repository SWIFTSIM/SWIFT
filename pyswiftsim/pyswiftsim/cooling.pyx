cimport pyswiftsim
from pyswiftsim cimport pointer
from pyswiftsim.pointer cimport Pointer
from pyswiftsim.pointer cimport pointer_check
from pyswiftsim.pointer cimport pointer_create

from pyswiftsim cimport part_def
from pyswiftsim.part_def cimport hydro_space

from pyswiftsim cimport equation_of_state
from pyswiftsim.equation_of_state cimport gas_entropy_from_internal_energy
from pyswiftsim.equation_of_state cimport gas_soundspeed_from_internal_energy

from libc.stdlib cimport malloc, free
from libc.math cimport M_PI
cimport numpy as np
import numpy as np

ctypedef np.float_t DTYPE_t

_hydro_gamma = hydro_gamma


def read_params(str filename):
    # transform python to C
    tmp = filename.encode(u"ascii")
    cdef char *cfilename = tmp
    if cfilename == NULL:
        raise Exception("Unable to convert filename to char*")

    # init params
    cdef swift_params *params = <swift_params*> malloc(sizeof(swift_params));
    parser_read_file(cfilename, params)
    #parser_print_params(params)
 
    
    return pointer_create(
        params,
        pointer.type_swift_params
    )

def init_units(Pointer p_params):

    # deal with inputs
    pointer_check(p_params, pointer.type_swift_params)
    cdef swift_params *params = <swift_params*> p_params._data
    
    # init units
    cdef unit_system *units = <unit_system*> malloc(sizeof(unit_system));
    cdef const char* category = "InternalUnitSystem"
    units_init(units, params, category)

    #units_print_backend(units)

    # init constants
    cdef phys_const *constants = <phys_const*> malloc(sizeof(phys_const));
    phys_const_init(units, constants)

    # return
    d = {}
    d["constants"] = pointer_create(constants, pointer.type_phys_const)
    d["units"] = pointer_create(units, pointer.type_unit_system)
    return d

def pycooling_init(Pointer p_params,
                   Pointer p_units,
                   Pointer p_constants):
    # deal with inputs
    pointer_check(p_params, pointer.type_swift_params)
    cdef swift_params *params = <swift_params*> p_params._data

    pointer_check(p_units, pointer.type_unit_system)
    cdef unit_system *units = <unit_system*> p_units._data

    pointer_check(p_constants, pointer.type_phys_const)
    cdef phys_const *constants = <phys_const*> p_constants._data

    #units_print_backend(units)
    # init cooling
    cdef cooling_function_data *cooling = <cooling_function_data*> malloc(sizeof(cooling_function_data));
    cooling_init(params, units, constants, cooling)

    # print results
    #cooling_print_backend(cooling)

    # store and return results
    return pointer_create(cooling, pointer.type_cooling_function_data)


def pycooling_rate(Pointer p_units,
                   Pointer p_cooling,
                   Pointer p_constants,
                   np.ndarray [DTYPE_t, ndim=1] rho,
                   np.ndarray[DTYPE_t, ndim=1] u):

    if rho.shape[0] != u.shape[0]:
        raise Exception("Rho and u should be of the same size!")
    
    cdef int N = rho.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] rate = np.empty_like(rho)
    
    # deal with inputs
    # p_units
    pointer_check(p_units, pointer.type_unit_system)
    cdef unit_system *units = <unit_system*> p_units._data

    # p_cooling
    pointer_check(p_cooling, pointer.type_cooling_function_data)
    cdef cooling_function_data *cooling = <cooling_function_data*> p_cooling._data

    # p_constants
    pointer_check(p_constants, pointer.type_phys_const)
    cdef phys_const *constants = <phys_const*> p_constants._data

    # initialize particles
    cdef part *p = <part*> malloc(sizeof(part) * N)
    cdef xpart *xp = <xpart*> malloc(sizeof(xpart) * N)

    cdef hydro_space hs # shadowfax stuff => no need to initialize

    for i in range(N):
        # init
        part_def.hydro_init_part(&p[i], &hs)
        p[i].id = i
        p[i].mass = 100.
        p[i].h = 1.
        p[i].rho = rho[i]
        p[i].entropy = gas_entropy_from_internal_energy(rho[i], u[i])
        for j in range(3):
            p[i].x[j] = 0
            p[i].v[j] = 0

    # compute cooling
    for i in range(N):
        # compute cooling
        rate[i] = cooling_rate(constants, units, cooling, &p[i])

        
    return rate


def soundspeed_from_internal_energy(np.ndarray [DTYPE_t, ndim=1] rho,
                                    np.ndarray [DTYPE_t, ndim=1] u):
    

    if rho.shape[0] != u.shape[0]:
        raise Exception("Rho and u should be of the same size!")
    
    cdef int N = rho.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] cs = np.empty_like(rho)

    for i in range(N):
        cs[i] = gas_soundspeed_from_internal_energy(rho[i], u[i])

    return cs
