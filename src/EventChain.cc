#include "Random.h"
#include "EventChain.h"
#include <cmath>

mcec::EventChain::EventChain() {}

void mcec::EventChain::evolveSystem(std::vector<mcec::Particle*>& particles, const unsigned int steps)
{
	py::print("Running: ", steps);

    int disk_displacement = 10;


    for (size_t i = 0; i < steps; ++i)
    {        
        mcec::Random::get()->setRandomPRNGSeed();
        int dirc = mcec::Random::get()->randBool();

        float distance_to_go = disk_displacement;

        Particle* next_particle = particles.at(mcec::Random::get()->randintRange(0, particles.size()-1));


        while(distance_to_go > 0.0) 
        {

            Particle* cpy_particle = next_particle;
            float min_particle_distance = distance_to_go;

            for (auto p : particles)
            {
                float event_b = event(p, cpy_particle, dirc);
                if (event_b < min_particle_distance)
                {
                    next_particle = p;
                    min_particle_distance = event_b;
                }
            }
            if (dirc == 1) next_particle->setCoordX(fmod((next_particle->getCoordX() + min_particle_distance), 1.0));
            else next_particle->setCoordY(fmod((next_particle->getCoordY() + min_particle_distance), 1.0));
            distance_to_go -= min_particle_distance;
        }

    }
}

float mcec::EventChain::event(mcec::Particle* p1, mcec::Particle* p2, int dirc)
{
    float d_perp = 0.0;
    float d_para = 0.0;

    if (dirc == 1)
        d_perp = std::fmod(fabs(p1->getCoordY() - p2->getCoordY()), 1.0);
    else
        d_perp = std::fmod(fabs(p1->getCoordX() - p2->getCoordX()), 1.0);

    d_perp = fmin(d_perp, 1.0 - d_perp);
    if (d_perp > 2.0 * p1->getRadius())
        return INFINITY;
    else
    {
        d_para = sqrt(fabs(4.0 * std::pow(p1->getRadius(), 2) - std::pow(d_perp, 2)));

        if (dirc == 1)
        {
            return std::fmod((p1->getCoordX() - p2->getCoordX() - d_para + 1.0), 1.0);
        }
        else
        {
            return std::fmod((p1->getCoordY() - p2->getCoordY() - d_para + 1.0), 1.0);
        }
    }
}