#include <random>
#include "EventChain.h"


mcec::EventChain::EventChain() {}

void mcec::EventChain::evolveSystem(std::vector<mcec::Particle*>& particles, unsigned int steps)
{
	py::print("Running: ", steps);

    int disk_displacement = 10;


    for (size_t i = 0; i < steps; ++i)
    {
        boost::random::uniform_int_distribution<> dist{1, 100};
        
        py::print("num = ", dist(m_gen));
    }
}

float mcec::EventChain::checkCollisionOfParticles(mcec::Particle p1, mcec::Particle p2)
{
    /*float d_perp = 0.0;
    float d_para = 0.0;

    if (dirc == 1)
        d_perp = std::fmod(fabs(b_y - a_y), 1.0);
    else
        d_perp = std::fmod(fabs(b_x - a_x), 1.0);

    d_perp = fmin(d_perp, 1.0 - d_perp);
    if (d_perp > 2.0 * sigma)
        return (float)INFINITY;
    else{
        d_para = sqrt(fabs(4.0 * sigma*sigma - d_perp*d_perp));

        if (dirc == 1){
            return std::fmod((b_x - a_x - d_para + 1.0), 1.0);
        }
        else{
            return std::fmod((b_y - a_y - d_para + 1.0), 1.0);
        }
    }*/
        return 23.2;
}