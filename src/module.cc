#include <pybind11/pybind11.h>
#include "System.h"
#include "Particle.h"

namespace py = pybind11;


PYBIND11_MODULE(mcec, m)
{

    init_system(m);
    init_particle(m);


        
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
    // return m.ptr();
}
