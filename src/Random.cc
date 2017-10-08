#include "Random.h"

mcec::Random* mcec::Random::m_instance = 0;

mcec::Random::Random()
{
    m_seed = 0;
    m_generator = std::mt19937(m_seed);

}

int16_t mcec::Random::randBool()
{
    std::uniform_int_distribution<> dis(0, 1);
    return dis(m_generator);
}

int mcec::Random::randintRange(int32_t a, int32_t b)
{
    std::uniform_int_distribution<> dis(a, b);
    return dis(m_generator);
}

float mcec::Random::random()
{
    std::uniform_real_distribution<> dis(0, 1);
    return dis(m_generator);
}
