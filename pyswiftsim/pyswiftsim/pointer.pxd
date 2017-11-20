cdef extern from "clocks.h":
    void clocks_set_cpufreq(unsigned long long freq);

    
cdef enum pointer_type:
    type_swift_params
    type_unit_system
    type_phys_const
    type_cooling_function_data
    type_count

cdef const char **pointer_name


cdef class Pointer:
    cdef void *_data
    cdef int _data_type
    cdef _setup(self, void *params, int data_type)


cdef pointer_create(void *params, int data_type)
cdef pointer_check(Pointer ptr, int data_type)
