cdef extern from "equation_of_state.h":
    cdef float gas_entropy_from_internal_energy(float density, float u)
    cdef float gas_soundspeed_from_entropy(float density, float entropy)
    cdef float gas_soundspeed_from_internal_energy(float density, float u)
