#include "System.h"
#include <cmath>

mcec::System::System(unsigned int num_particles, float density, float radius) 
{
	m_num_particles = num_particles;
	m_density = density;
	m_radius = radius;
	m_disk_displacement = std::sqrt(m_num_particles); 
}

void init_system(pybind11::module& m)
{
	pybind11::class_<mcec::System>(m, "System")
        .def(pybind11::init<unsigned int, float, float>(),py::arg("num_particles")=0, py::arg("density")=0., py::arg("radius")=0.)
		.def("setNumParticles", &mcec::System::setNumParticles)
		.def("getNumParticles", &mcec::System::getNumParticles)
		.def("setDensity", &mcec::System::setDensity)
		.def("getDensity", &mcec::System::getDensity)
		.def("getParticleRadius", &mcec::System::getParticleRadius);
}

