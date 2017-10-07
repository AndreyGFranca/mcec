#pragma once
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>

#include "Particle.h"
#include "EventChain.h"

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
		// Eigen::MatrixXd getParticleList() { return m_particle_positions; }
		std::vector<Particle*> getParticleList() { return m_particles; }

		void init();
		void run(unsigned int steps);



	private:
		unsigned int m_num_particles;
		float m_density;
		float m_radius;
		float m_disk_displacement;
		float m_delxy;
		float m_two_delxy;
		std::string m_initial_configuration;
		// Eigen::MatrixXd m_particle_positions;
		std::vector<Particle*> m_particles;
		std::shared_ptr<mcec::EventChain> m_event_chain;

	};

}

void init_system(pybind11::module& m);


