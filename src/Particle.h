#pragma once

#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace mcec
{

	class Particle
	{
	public:
		Particle(float radius, float x, float y);

		float getRadius() const { return m_radius; }
		float getCoordX() const { return m_x; }
		float getCoordY() const { return m_y; }
		void setCoordX(float x) { m_x = x; }
		void setCoordY(float y) { m_y = y; }

		bool operator == (const Particle& rhs)
		{
			return ( m_x == rhs.getCoordX() && m_y == rhs.getCoordY() );
		}
		// ~Particle();

	private:
		float m_radius;
		float m_x, m_y;
		float m_density;

	};

}

void init_particle(py::module& m);