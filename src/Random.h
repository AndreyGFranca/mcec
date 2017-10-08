#include <random>
#include <cstdint>
#include <chrono>

namespace mcec
{
    class Random
    {

    public:
        static Random* get()
        {
            if(m_instance == 0)
            {
                m_instance = new Random();
                return m_instance;
            }
            return m_instance;
        }
        ~Random();

        int64_t setPRNGSeed(int64_t seed)
        {
            m_seed = seed;
            m_generator.seed(seed);
        }

        void setRandomPRNGSeed()
        {
            std::random_device r;
            std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
            m_generator.seed(seed);
        }

        int16_t randBool();
        float random();

        int randintRange(int32_t a, int32_t b);



    private:
        Random();
        std::random_device rd;
        static Random* m_instance;
        int64_t m_seed;
        std::mt19937 m_generator;        
    };

}