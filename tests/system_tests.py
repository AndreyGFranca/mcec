import mcec
# from mcec import System

t = mcec.System(10, 5.0)
assert ( t.getNumParticles() == 10 )

t.setNumParticles(20)

assert(t.getNumParticles() == 20)

t.setDensity(20.32)
print(t.getDensity())


t2 = mcec.System(num_particles=8**2, density=5.0, radius=2)

a = t2.getParticleList()

print (a)

t2.init()
a = t2.getParticleList()

print (a)





print(t)
