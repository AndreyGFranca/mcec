#pragma once
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

namespace mcec
{
	class System
	{
	public:
		System(unsigned int num_particles, float density, float radius);
		//~System();

		void setNumParticles(unsigned int num_particles) { m_num_particles = num_particles; }
		unsigned int getNumParticles() const { return m_num_particles; }
    	void setDensity(float density) { m_density = density; }
		float getDensity() const { return m_density; }
		float getParticleRadius() const { return m_radius; }


	private:
		unsigned int m_num_particles;
		float m_density;
		float m_radius;
		float m_disk_displacement;
		std::string m_initial_configuration;

	};

}

void init_system(pybind11::module& m);


