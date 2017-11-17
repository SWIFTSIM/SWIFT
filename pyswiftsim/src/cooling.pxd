
cdef extern from "clocks.h":
    void clocks_set_cpufreq(unsigned long long freq);
    
cdef extern from "parser.h":
    struct swift_params:
        pass

    void parser_read_file(const char *file_name, swift_params *params);
    void parser_print_params(const swift_params *params);

cdef extern from "units.h":
    struct unit_system:
        pass
    void units_init(unit_system* us, const swift_params* params,
                    const char* category)

    
cdef extern from "physical_constants.h":
    struct phys_const:
        pass
    void phys_const_init(unit_system* us, phys_const* internal_const);

    
cdef extern from "cooling_struct.h":
    struct cooling_function_data:
        pass

cdef extern from "part.h":
    struct part:
        pass
    struct xpart:
        pass

cdef extern from "hydro_space.h":
    struct hydro_space:
        pass

cdef extern from "hydro.h":
    void hydro_init_part(
        part *p, const hydro_space *hs)

cdef extern from "debug.h":
    void printParticle(const part *parts, const xpart *xparts,
                       long long int id, size_t N);

cdef extern from "cooling.h":
    void cooling_init(const swift_params* parameter_file,
                      const unit_system* us,
                      const phys_const* phys_const,
                      cooling_function_data* cooling);
    void cooling_print_backend(
        const cooling_function_data* cooling)

    void cooling_cool_part(
        const phys_const* phys_const,
        const unit_system* us,
        const cooling_function_data* cooling,
        part* p, xpart* xp, float dt)
