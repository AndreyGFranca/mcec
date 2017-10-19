import numpy as np
import matplotlib.pyplot as plt
import math, time
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#####################################
###Funcoes

def exponencial(x, A0, A1, tau):
	f = A0+A1*np.exp(-x/tau)
	return f

def ajuste_exp(eixox, eixoy, chute):
	parametros, erros = curve_fit(exponencial, eixox, eixoy, chute)
	erroA0 = np.sqrt(erros[0,0])
	erroA1 = np.sqrt(erros[1,1])
	errotau = np.sqrt(erros[2,2])
	return parametros, [erroA0, erroA1, errotau]

def grafico_ajusteC(tempo,cdelta,ajuste,nome,ind,tau_erro):
	plt.plot(tempo, cdelta, 'ob', label = r'$C(\delta)$')
	plt.plot(tempo, ajuste, '-r', label = r'$e^{-\delta/\tau}$')
	titulo = r'Exponential fit: $\eta=$ '+str(eta)+', $p=$ '+str(lista_n0[ind])
	plt.title(titulo,fontsize=18)
	plt.xticks(color='k', size=12)
	plt.yticks(color='k', size=12)
	plt.legend(loc='upper right',fontsize = 16)
	plt.xlabel(r'Time $\delta$', fontsize = 16) #k indica cor preta
	plt.ylabel('Correlation function', fontsize = 16)
	plt.grid(True)
	#axes = plt.gca()
	#axes.set_xlim([5,10])
	#axes.set_ylim([-1.0,-0.7])
	plt.text(15000,0.7,tau_erro,fontsize = 14, backgroundcolor = 'w', color = 'k')
	#plt.savefig(nome, dpi = 300)
	plt.close()

############################################3
###Importando os dados

#Inicio da marcacao do tempo de calculo
start_time = time.clock()

Q = 25000
samples = 5
per = 100
eta = 0.705
steps = int(Q/per)

lista_n0 = [64, 128, 256, 512]
Nn0 = len(lista_n0)

tempo = np.empty(steps)
Psi = np.empty((Nn0,steps))
Psi2 = np.empty((Nn0,steps))

i1 = 0
for n0 in lista_n0:
	fname = 'Psi_hex_n0'+str(n0)+'_eta'+str(int(1000*eta))+'_t'+str(Q)+'_S'+str(samples)+'_p'+str(per)+'.csv'
	file = open(fname,"r")
	t1 = 0
	for linha in file:
		a1, b1, c1 = linha.split(',')
		tempo[t1] = int(a1)
		Psi[i1][t1] = float(b1) / samples
		Psi2[i1][t1] = float(c1) / samples  
		t1 += 1
	file.close()
	i1 += 1


Corr = np.zeros((Nn0,steps))

for i3 in range(Nn0):
	Psi[i3] = Psi[i3] - Psi[i3][steps-1] 
	Psi2[i3] = Psi2[i3] - Psi2[i3][steps-1] 
	for delta in range(steps):
		for t3 in np.arange(0,steps-1-delta,1):
			Corr[i3][delta] += Psi[i3][t3]*Psi[i3][t3+delta] 
	Corr[i3] = Corr[i3]/Corr[i3][0] #normalizando

##########################################

#####################################################
#Graficos

titulo2 = r'Correlation Function: $\eta=$ '+str(eta)+', $S=$ '+str(samples)+'.'

for i2 in range(Nn0):
	plt.plot(tempo, Corr[i2], cores[i2], lw = 2, label = r'$N=$ '+str(lista_n0[i2])+'$^2$')
#plt.plot(tempo, ajuste, '-g', label = r'$f(x)$')
#plt.title(titulo,fontsize=16)
plt.xticks(color='k', size=12)
plt.yticks(color='k', size=12)
plt.legend(loc='upper right',fontsize = 16)
plt.xlabel(r'Time $\delta$', fontsize = 16) #k indica cor preta
plt.ylabel('$C(\delta)$', fontsize = 16)
plt.grid(True)
#axes = plt.gca()
#axes.set_xlim([5,10])
#axes.set_ylim([-1.0,-0.7])
plt.text(17000,0.5,r'$\eta=$ '+str(eta),fontsize = 16, backgroundcolor = 'w', color = 'k')
plt.savefig('figuraCorr.pdf', dpi = 300)
plt.close()

##########Ajuste

chuteC = []

chuteC0 = [-7.4554e-02, 1.023, 3.795e+03] #initial values
chuteC.append(chuteC0)

chuteC1 = [-8.79e-03,8.765e-01, 6.96e+03] #initial values
chuteC.append(chuteC1)

chuteC2 = [-4.982e-02, 9.981e-01, 7.562e+03] #initial values
chuteC.append(chuteC2)

chuteC3 = [-4.982e-02, 9.981e-01, 7.562e+03] #initial values
chuteC.append(chuteC3)

paraC, errosC, ajusteC = [], [], []

for i5 in range(Nn0):
	parametros, erros = ajuste_exp(tempo, Corr[i5], chuteC[i5])
	paraC.append(parametros)
	errosC.append(erros)
	ajusteC.append(paraC[i5][0] + paraC[i5][1]*np.exp(-tempo/paraC[i5][2]))
	print('p =',lista_n0[i5])
	print('Parametros',paraC[i5])
	print('Erros',errosC[i5])
	nome = 'ajuste_Corr_p'+str(round(lista_n0[i5]))+'.png'
	tau_erro = r'$\tau=$ '+str(int(paraC[i5][2]))+r' $\pm$ '+str(int(errosC[i5][2]))
	grafico_ajusteC(tempo,Corr[i5],ajusteC[i5],nome,i5,tau_erro)


filename2 = 'tau_Corr_eta'+str(int(1000*eta))+'_S'+str(samples)+'_p'+str(per)+'.csv'
arq2 = open(filename2, 'w')
for i6 in range(Nn0):
	arq2.write(str(lista_n0[i6])+','+str(paraC[i6][0])+','+str(paraC[i6][1])+','+str(round(paraC[i6][2]))+','+str(round(errosC[i6][2],2))+'\n')
arq2.close() 

