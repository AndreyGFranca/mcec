#include "EventChain.h"

mcec::EventChain::EventChain() {}

void mcec::EventChain::evolveSystem(std::vector<mcec::Particle*>& particles, unsigned int steps)
{
	py::print("Running: ", steps);

    int disk_displacement = 10;

    for (size_t i = 0; i < steps; ++i)
    {
        

    }
}