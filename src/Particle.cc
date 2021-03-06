#include "Particle.h"
#include <string>


mcec::Particle::Particle(float radius, float x, float y) 
{
	m_radius = radius;
	m_x = x;
	m_y = y;
}

void init_particle(py::module& m)
{
	py::class_<mcec::Particle>(m, "Particle")
		.def(py::init<float, float, float>(),py::arg("radius")=0., py::arg("x")=0., py::arg("y")=0.)
		.def("getRadius", &mcec::Particle::getRadius)
		.def_property("x", &mcec::Particle::getCoordX, &mcec::Particle::setCoordX)
		.def_property("y", &mcec::Particle::getCoordY, &mcec::Particle::setCoordY)
        .def("__getitem__", [](const mcec::Particle &p) {
                return std::pair<float, float> (p.getCoordX(), p.getCoordY());
            })
        .def("__str__",
            [](const mcec::Particle &p) {
                return "(" + std::to_string(p.getCoordX()) + ", " + std::to_string(p.getCoordY())+ ")";
            }
        )
        .def("__repr__",
            [](const mcec::Particle &p) {
                return "(" + std::to_string(p.getCoordX()) + ", " + std::to_string(p.getCoordY())+ ")";
            }
        );
}
