import hoomd
import hoomd.hpmc
import freud
import numpy as np
from freud.order import HexOrderParameter
import matplotlib.pyplot as plt
import math, time#, random

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class PsiAnalyzer:

	def __init__(self, system):
		self.system = system
		self.psi_soma_array = np.zeros(Q, dtype=np.complex64)
		self.psi2_soma_array = np.zeros(Q, dtype=np.complex64)
		self.order_parameter = HexOrderParameter(2.8, 6.0, 6)
		#self.total_particles = system.take_snapshot().particles.N
		self.initial_configuration = system.take_snapshot();

	def __call__(self, step):
		snap = self.system.take_snapshot()
		pos = snap.particles.position
		box = freud.box.Box(snap.box.Lx, snap.box.Ly, is2D=True)
		self.order_parameter.compute(box, pos)
		self.psi_soma_array[step % Q] += np.sum(self.order_parameter.getPsi())
		self.psi2_soma_array[step % Q] += np.sum(self.order_parameter.getPsi()**2)


	#Metodos
	def get_soma_psi_abs(self):
		return  np.abs(self.psi_soma_array)

	def get_soma_psi2_abs(self):
		return  np.abs(self.psi2_soma_array)

	def get_soma_psi_cplx(self):
		return  self.psi_soma_array

	def get_soma_psi_real(self):
		return  self.psi_soma_array.real

	def get_soma_psi_imag(self):
		return  self.psi_soma_array.imag

start_time = time.clock() #Time count

#lista_eta = [705, 708, 711, 714, 717, 720, 723, 735]

#Parametros
n0 = 256
L = 100.0
twodelxy = L/n0 #distancia entre dois discos consecutivos 
eta = 0.705
Q = 25000
samples = 5
per = 100 #periodicidade, calcula psi a cada per iteracoes

vetor_Psi, vetor_Psi2 = np.zeros(Q), np.zeros(Q)

for i in range(samples):
	hoomd.context.initialize("--mode=cpu")
	system = hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=twodelxy), n=n0);
	#print(system.box)
	area = system.box.Lx*system.box.Ly
	snapshot = system.take_snapshot(all=True)
	analyzer = PsiAnalyzer(system); #criando uma instancia
	hoomd.analyze.callback(analyzer, period = per)
	Np = snapshot.particles.N #numero de particulas
	radius = math.sqrt((area*eta) / (math.pi * Np))
	mc = hoomd.hpmc.integrate.sphere(d=0.2, seed=i) #definicao do formato do disco
	mc.shape_param.set('A', diameter=2*radius); #diametro do disco
	hoomd.run(Q,quiet = False) #evolucao temporal do sistem em Q iteracoes
	#calculando Psi
	analyzer.psi_soma_array = analyzer.psi_soma_array / Np
	#calculando Psi ao quadrado
	analyzer.psi2_soma_array = analyzer.psi2_soma_array / Np
	vetor_Psi += analyzer.get_soma_psi_abs()  #vetor com os valores de Psi global para cada instante de tempo
	vetor_Psi2 += analyzer.get_soma_psi2_abs()  #vetor com os valores de Psi global para cada instante de tempo


######## Gravacao dos resultados em um arquivo externo ##############

tempo = np.arange(Q)
filenome2 = 'Psi_hex_n0'+str(n0)+'_eta'+str(int(1000*eta))+'_t'+str(Q)+'_S'+str(samples)+'_p'+str(per)+'.csv'
arquivo2 = open(filenome2, 'w') #exportar a configuracao inicial
for i2 in np.arange(0,Q,per):
	arquivo2.write(str(tempo[i2])+' , '+str(vetor_Psi[i2])+' , '+str(vetor_Psi2[i2])+'\n')
arquivo2.close() 


tempototal = round(time.clock() - start_time,4)
tempo_medio = round(tempototal/ (Q*samples),4)
print("Average time per iteration: ",tempo_medio, "seconds")
print("Total execution time: ",tempototal, "seconds")

filenome2 = 'parametros_n0'+str(n0)+'_eta'+str(int(1000*eta))+'_t'+str(Q)+'_S'+str(samples)+'_p'+str(per)+'.txt'
arquivo2 = open(filenome2, 'w') 
arquivo2.write('Parametros. \n')
arquivo2.write('n0: '+str(n0)+'\n')
arquivo2.write('Numero de discos: '+str(Np)+'\n')
arquivo2.write('Lados da caixa Lx e Ly: '+str(system.box.Lx)+' ,'+str(round(system.box.Ly,3))+'\n')
arquivo2.write('Raio do disco: '+str(radius)+'\n')
arquivo2.write('Densidade eta: '+str(eta)+'\n')
arquivo2.write('Numero de iteracoes: '+str(Q)+'\n')
arquivo2.write('Numero de amostras: '+str(samples)+'\n')
arquivo2.write('Periodicidade: '+str(per)+'\n')
arquivo2.write('Modo de calculo: cpu \n')
arquivo2.write('Computador: phonon (Zmax) \n')
arquivo2.write('Tempo total de calculo em segundos: '+str(tempototal)+'\n')
arquivo2.write('Tempo total de calculo em horas: '+str(round(tempototal/3600.0,3))+'\n')
arquivo2.write('Tempo calculo por iteracao em segundos: '+str(tempo_medio)+'\n')
arquivo2.close()

print("The End.")



#############################################################################






