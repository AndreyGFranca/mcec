#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <vector>

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

    };

}