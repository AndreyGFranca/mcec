from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy import spatial
import numpy as np
import matplotlib.pyplot as plt
import math, random, pylab, time, cmath, matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def delx_dely(x, y):
	d_x = (x[0] - y[0]) % 1.0
	if d_x > 0.5: d_x -= 1.0
	d_y = (x[1] - y[1]) % 1.0
	if d_y > 0.5: d_y -= 1.0
	return d_x, d_y

def dist(x,y): #x e y sao lista ou vetores com duas coordenadas cada um
	d_x = abs(x[0] - y[0]) % 1.0
	d_x = min(d_x, 1.0 - d_x)
	d_y = abs(x[1] - y[1]) % 1.0
	d_y = min(d_y, 1.0 - d_y)
	return  math.sqrt(d_x**2 + d_y**2)


Nsqrt = 128
Q = 50
eta = 0.7
instante = Q

N = Nsqrt ** 2 #numero de discos, tem que ser um numero cuja raiz quadrada seja inteira
sigma = math.sqrt(eta / (N * math.pi)) #raio do disco em funcao de N e da densidade eta

posfinal = open("pos_final_128.csv", "r")
posicoes = []
for linha in posfinal:
	aQ, bQ = linha.split(',')
	posicoes.append([float(aQ), float(bQ)])

###################################################################

Nall = len(posicoes)
tree = spatial.KDTree(posicoes)

Psi = 0.0j
psik_real = np.empty(Nall)
psik_imag = np.empty(Nall)

for i in range(Nall):
	vector, n_neighbor  = 0.0j, 0
	discoi = np.array(posicoes[i]) #defino o disco para se calcular psi
	vizinhos = tree.query(discoi,6) #determina os 6 vizinhos mais proximos
	vizinhosID = np.array(vizinhos[1]) #vetor com os indices dos 6 vizinhos mais proximos
	for j in vizinhosID:
		dx, dy = delx_dely(posicoes[j], posicoes[i])
		angle = cmath.phase(complex(dx, dy))
		vector += cmath.exp(6.0j * angle) 
		n_neighbor += 1
	if n_neighbor > 0:
		vector /= float(n_neighbor)
	Psi += vector / float(Nall)
	psik_real[i] = vector.real 
	psik_imag[i] = vector.imag 

Psi_mod = abs(Psi) 
Psi_real = Psi.real
Psi_imag = Psi.imag

costheta = np.empty(Nall)
norma = np.empty(Nall)

#definicao do cosseno do angulo em funcao do disco
for k in range(Nall):
    norma[k] = Psi_mod * math.sqrt( psik_real[k] ** 2 + psik_imag[k] ** 2 )
    costheta[k] = (psik_real[k] * Psi_real + psik_imag[k] * Psi_imag) / norma[k]
    costheta[k] = round(costheta[k],4)

# ######################################################################
#Diagrama de Voronoi
vor = Voronoi(np.array(posicoes))
voronoi_plot_2d(vor)
vertices = []
nr = -1
for region in vor.regions:
	if not -1 in region:
		nr += 1
		poligono = [vor.vertices[i] for i in region] #selecionando os vertices de cada regiao
		vertices.append(poligono) #juntando os vertices de cada regiao em um poligono
		# plt.fill(*zip(*polygon))

if [] in vertices: #removendo as listas vazias
	vertices.remove([])

N = len(vertices)

# ########################################################################33

fig, ax = plt.subplots()
patches = []

for i in range(N):
	polygon = Polygon(vertices[i], True)
	patches.append(polygon) #fazendo os vertices um poligono

p = PatchCollection(patches, cmap=matplotlib.cm.magma, alpha=1.0, lw = 0.0)

p.set_array(costheta)
ax.add_collection(p)

# #######################################################################

#Grafico

Psim = round(Psi_mod,2)

nome3 = 'costhetak'+'_N'+str(Nsqrt)+'sq_t'+str(instante) + '_eta'+str(eta) + '.png'
tit = r'$N=$ '+str(Nsqrt)+'$^2$,  $\eta = $ '+str(eta)+',  $t = $ '+str(instante)+r',  $\vert \Psi \vert =$ '+str(Psim)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
clb = plt.colorbar(p)
clb.set_label(r'$\cos \theta_k$', fontsize = 18, labelpad=-40, y=1.05, rotation=0)
plt.axis('scaled')
plt.xlabel('eixo x')
plt.ylabel('eixo y')
plt.title(tit)
plt.axis([0.0, 1.0, 0.0, 1.0])
#plt.savefig(nome3,dpi=200)
plt.savefig('voronoi_128.pdf',dpi = 300)
#plt.show()

