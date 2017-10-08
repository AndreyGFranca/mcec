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

        void evolveSystem(std::vector<mcec::Particle*>& particles, const unsigned int steps);
        float event(mcec::Particle* p1, mcec::Particle* p2, int dirc);

    private:
        int seed;
        
    };

}