#pragma once
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <string>
#include <vector>

#include "Particle.h"

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
		Eigen::MatrixXd getParticleList() { return m_particle_positions; }

		void init();


	private:
		unsigned int m_num_particles;
		float m_density;
		float m_radius;
		float m_disk_displacement;
		float m_delxy;
		float m_two_delxy;
		std::string m_initial_configuration;
		Eigen::MatrixXd m_particle_positions;

	};

}

void init_system(pybind11::module& m);


