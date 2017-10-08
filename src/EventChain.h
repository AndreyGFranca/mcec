#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <vector>
#include <cstdint>
#include <random>

#include "Particle.h"

namespace py = pybind11;

namespace mcec
{

	class EventChain
    {
    public:
        EventChain();
        ~EventChain();

        void evolveSystem(std::vector<mcec::Particle*>& particles, unsigned int steps);
        float checkCollisionOfParticles(mcec::Particle p1, mcec::Particle p2);

    private:
        int seed;
        boost::random::mt19937 m_gen{static_cast<std::uint32_t>(1)};
    };

}