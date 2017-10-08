import mcec 

t = mcec.System(num_particles=8**2, density=4.0)

assert(t.getNumParticles() == 64)

print (t.getParticleList())

t.init()

#print (t.getParticleList())

t.run(20)
print (t.getParticleList())