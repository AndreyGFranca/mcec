import mcec
# from mcec import System

t = mcec.System(10, 5.0)
assert ( t.getNumParticles() == 10 )

t.setNumParticles(20)

assert(t.getNumParticles() == 20)

t.setDensity(20.32)
print(t.getDensity())


t2 = mcec.System(num_particles=10, density=5.0, radius=2)

assert(t2.getParticleRadius() == 2)


print(t)
