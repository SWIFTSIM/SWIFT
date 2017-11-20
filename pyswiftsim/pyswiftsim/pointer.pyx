from libc.stdlib cimport free
from cooling cimport unit_system, phys_const

clocks_set_cpufreq(0); # estimate automatically cpufreq and init time

pointer_name = [
    "swift_params",
    "unit_system",
    "phys_const",
    "cooling_function_data",
    "count"
]

# wrapper around C pointers
cdef class Pointer:
    def __cinit__(self):
        self._data = NULL

    cdef _setup(self, void *params, int data_type):
        self._data = params
        self._data_type = data_type
        return self

    def __str__(self):
        txt = "%s" % pointer_name[self._data_type]
        return txt

    def __del__(self):
        free(self._data)

    def __getattr__(self, value):
        data_type = self._data_type
        if data_type == type_unit_system:
            return pointer_unit_system_getattr(self, value)
        if data_type == type_phys_const:
            return pointer_phys_const_getattr(self, value)
        else:
            raise Exception("Data type not implemented")


# method to call when creating a pointer
cdef pointer_create(void *params, int data_type):
    return Pointer()._setup(params, data_type)

# method to call when receiving a pointer
cdef pointer_check(Pointer ptr, int data_type):
    if (ptr._data_type != data_type):
        exp = pointer_name[data_type]
        rec = pointer_name[ptr._data_type]
        raise Exception("Expecting pointer of type %s, got %s" % (exp, rec))
    if (ptr._data == NULL):
        data_name = pointer_name[ptr._data_type]
        raise Exception("Value not set for pointer type %s" % data_name)
    

cdef pointer_unit_system_getattr(Pointer ptr, str attr):
    pointer_check(ptr, type_unit_system)
    cdef unit_system us = (<unit_system*> ptr._data)[0]

    if attr == "UnitMass_in_cgs":
        return us.UnitMass_in_cgs
    elif attr == "UnitLength_in_cgs":
        return us.UnitLength_in_cgs
    elif attr == "UnitTime_in_cgs":
        return us.UnitTime_in_cgs
    elif attr == "UnitCurrent_in_cgs":
        return us.UnitCurrent_in_cgs
    elif attr == "UnitTemperature_in_cgs":
        return us.UnitTemperature_in_cgs
    else:
        raise Exception("unit_system does not contain %s" % attr)

cdef pointer_phys_const_getattr(Pointer ptr, str attr):
    pointer_check(ptr, type_phys_const)
    cdef phys_const consts = (<phys_const*> ptr._data)[0]

    if attr == "const_boltzmann_k":
        return consts.const_boltzmann_k
    if attr == "const_proton_mass":
        return consts.const_proton_mass
    else:
        raise Exception("phys_const does not contain %s" % attr)
