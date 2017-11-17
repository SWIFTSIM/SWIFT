cdef extern from "equation_of_state.h":
    cdef float gas_entropy_from_internal_energy(float density, float u)
