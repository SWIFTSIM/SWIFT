cdef extern from "part.h":
    struct part:
        float h
        float mass
        float rho
        long long id
        double x[3]
        float v[3]
        float entropy
        
    struct xpart:
        pass

cdef extern from "hydro.h":
    void hydro_init_part(
        part *p, const hydro_space *hs)

cdef extern from "hydro_space.h":
    struct hydro_space:
        pass

cdef extern from "debug.h":
    void printParticle(const part *parts, const xpart *xparts,
                       long long int id, size_t N);

