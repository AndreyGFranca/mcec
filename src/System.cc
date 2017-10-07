#include "System.h"
#include <cmath>

mcec::System::System(unsigned int num_particles, float density, float radius) 
{
	m_num_particles = num_particles;
	m_density = density;
	m_radius = std::sqrt(m_density /(m_num_particles* M_PI));
	m_disk_displacement = std::sqrt(m_num_particles); 
	m_initial_configuration = "quad";
	m_delxy = 1.0/ (2.0 * std::sqrt(m_num_particles));
	m_two_delxy = 2.0 * m_delxy;
	m_particle_positions = Eigen::MatrixXd::Zero(std::sqrt(m_num_particles), std::sqrt(m_num_particles));

}

void mcec::System::init()
{
	for (size_t i = 0; i < std::sqrt(m_num_particles); ++i)
	{
		for (size_t j = 0; j < std::sqrt(m_num_particles); ++j)
		{
			m_particle_positions(i, j) = m_delxy + j * m_two_delxy;
		}
	}
}

void init_system(pybind11::module& m)
{
	pybind11::class_<mcec::System>(m, "System")
        .def(pybind11::init<unsigned int, float, float>(),py::arg("num_particles")=0, py::arg("density")=0., py::arg("radius")=0.)
		.def("setNumParticles", &mcec::System::setNumParticles)
		.def("getNumParticles", &mcec::System::getNumParticles)
		.def("setDensity", &mcec::System::setDensity)
		.def("getDensity", &mcec::System::getDensity)
		.def("getParticleRadius", &mcec::System::getParticleRadius)
		.def("getParticleList", &mcec::System::getParticleList)
		.def("init", &mcec::System::init);
}

