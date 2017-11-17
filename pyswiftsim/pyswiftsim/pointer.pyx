from libc.stdlib cimport free

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
    

