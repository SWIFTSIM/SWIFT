clocks_set_cpufreq(0); # estimate automatically cpufreq and init time

# pointer types
cdef str _SWIFT_PARAMS = "swift_params"
cdef str _UNIT_SYSTEM = "unit_system"
cdef str _PHYS_CONST = "phys_const"
cdef str _COOLING_FUNCTION_DATA = "cooling_function_data"

# wrapper around C pointers
cdef class Pointer:
    cdef void *_data
    cdef str _data_type

    def __cinit__(self):
        self._data = NULL

    cdef _setup(self, void *params, str data_type):
        self._data = params
        self._data_type = data_type
        return self

    def __str__(self):
        return self._data_type


cdef pointer_create(void *params, str data_type):
    return Pointer()._setup(params, data_type)

cdef check_pointer(Pointer ptr, str data_type):
    if (ptr._data_type != data_type):
        raise Exception("Expecting pointer of type %s, got %s" % (data_type, ptr._data_type))
    if (ptr._data == NULL):
        raise Exception("Value not set for pointer type %s" % ptr._data_type)
    

def read_params(str filename):
    tmp = filename.encode(u"ascii")
    cdef char *cfilename = tmp
    if cfilename == NULL:
        raise Exception("Unable to convert filename to char*")

    # init params
    cdef swift_params params;
    parser_read_file(cfilename, &params)
    parser_print_params(&params)

    
    return pointer_create(&params, _SWIFT_PARAMS)

def init_units(Pointer p_params):

    # deal with inputs
    check_pointer(p_params, _SWIFT_PARAMS)
    cdef swift_params *params = <swift_params*> p_params._data
    
    # init units
    cdef unit_system units;
    cdef const char* category = "InternalUnitSystem"
    units_init(&units, params, category)

    # init constants
    cdef phys_const constants;
    phys_const_init(&units, &constants)

    # return
    d = {}
    d["constants"] = pointer_create(&constants, _PHYS_CONST)
    d["units"] = pointer_create(&units, _UNIT_SYSTEM)
    return d

def pycooling_init(Pointer p_params,
                   Pointer p_units,
                   Pointer p_constants):
    # deal with inputs
    check_pointer(p_params, _SWIFT_PARAMS)
    cdef swift_params *params = <swift_params*> p_params._data

    check_pointer(p_units, _UNIT_SYSTEM)
    cdef unit_system *units = <unit_system*> p_units._data

    check_pointer(p_constants, _PHYS_CONST)
    cdef phys_const *constants = <phys_const*> p_constants._data

    # init cooling
    cdef cooling_function_data cooling;
    cooling_init(params, units, constants, &cooling)

    # print results
    cooling_print_backend(&cooling)

    # store and return results
    return pointer_create(&cooling, _COOLING_FUNCTION_DATA)


def pycooling_rate(Pointer p_params):
    # deal with inputs
    check_pointer(p_params, _SWIFT_PARAMS)
    cdef swift_params *params = <swift_params*> p_params._data

    cdef part p
    cdef xpart xp

    cdef hydro_space hs # shadowfax stuff => no need to initialize
    hydro_init_part(&p, &hs)

    printParticle(&p, &xp, 0, 1)
    
    
