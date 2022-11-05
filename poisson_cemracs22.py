import numpy as np
import matplotlib.pyplot as plt

#domain is [-1,1]

#we solve partial_xE = rho
#with E 

def rho_fun0(x):
	return np.sin(np.pi*x)

def rho_fun1(x):
	return np.cos(np.pi*x)

N=10

eps = 2./N
#a = 0.1

a = -1.+3.*eps

def rho_fun2(x):
	if(np.abs(x)>eps):
		return 0.
	else:
		return (1.-np.abs(x)/eps)/eps	

def rho_fun3(x):
	if(np.abs(x-a)>eps):
		return 0.
	else:
		return (1.-np.abs(x-a)/eps)/eps	


def E_fun0(x):
	return (1.-np.cos(np.pi*x))/np.pi

def E_fun1(x):
	return np.sin(np.pi*x)/np.pi
	
def E_fun2(x):
	if(x>0):
		return 0.5
	else:
		if(x==0):
			return 0.
		else:
			return -0.5		

def E_fun3(x):
	if(x>a):
		return 0.5
	else:
		if(x==a):
			return 0.
		else:
			return -0.5		
		

def rho_fun(x):
	return rho_fun3(x)

def E_fun(x):
	return E_fun3(x)


def compute_rho(x,rho_fun):
	N = len(x)-1
	h = 2.0/N
	rho = np.zeros(N+1)
	for i in range(N+1):
		rho[i] = rho_fun(x[i])
	return rho

#we solve partial_xE = rho
#with E(0) = 0 
#     int i=0, i0 = (int)(N/2); // middle of the domain 
#     double factor = 0.5 * dx / (lambda*lambda); // 0.5 from trapezes
# 
#     if (2*i0 != N) {
#         printf("\n\n[update_E_from_rho_and_current_1d] WARNING : 0 n'appartient pas au maillage. Milieu i0 = %d, N=%d.\n\n", i0, N);
#     }
#     E[i0] = 0.0; // not even needed since it is never updated and E[0]=0 at t=0, but more safe
#     for (i=i0+1; i<=N; ++i) {
#         E[i] = E[i-1] + factor * (rho[i-1] + rho[i]);
#     }
#     for (i=i0-1; i>=0; --i) {
#         E[i] = E[i+1] - factor * (rho[i+1] + rho[i]);
#     }


def compute_E(rho):
	N = len(x)-1
	h = 2.0/N
	E = np.zeros(N+1)
	i0 = N//2
	factor = 0.5*h
	for  i in range(i0+1,N+1):
		E[i] = E[i-1]+factor*(rho[i-1]+rho[i])
	for  i in range(i0-1,-1,-1):
		E[i] = E[i+1]-factor*(rho[i+1]+rho[i])
	return E	

def compute_E2(rho):
	N = len(x)-1
	h = 2.0/N
	E = np.zeros(N+1)
	factor = 0.5*h
	C = 0.
	for  i in range(1,N+1):
		val = factor*(rho[i-1]+rho[i])
		E[i] = E[i-1]+val
		C += val
	print(C)	
	C = -C/2.
	for i in range(N+1):
		E[i] +=C	
	return E	



def compute_fun(x,fun):
	N = len(x)-1
	f = np.zeros(N+1)
	for i in range(N+1):
		f[i] = fun(x[i])
	return f
	
def compute_norm(tab,x):
	N = len(x)-1
	h = 2.0/N
	linf = np.max(np.abs(tab))
	val = np.abs(tab[0])
	l1 = 0.5*val
	l2 = 0.5*val*val 
	for i in range(1,N):
		val = np.abs(tab[i])
		l1 += val
		l2 += val*val 
	val = np.abs(tab[N])
	l1 += 0.5*val
	l2 += 0.5*val*val 
	l2 = np.sqrt(h*l2)
	l1 = h*l1
	return l1,l2,linf
	
x = np.linspace(-1,1,N+1)
rho = compute_rho(x,rho_fun)
E =  compute_E(rho)
E2 =  compute_E2(rho)
E_exact  = compute_fun(x,E_fun)

plt.plot(x,rho)
plt.show()

plt.plot(x,E,label='E(0)=0')
plt.plot(x,E2,label='int E = 0')
plt.plot(x,E_exact,label='E exact')
plt.legend()
plt.show()

err = E_exact-E
l1,l2,linf = compute_norm(err,x)
print(N,l1,l2,linf)
err = E_exact-E2
l1,l2,linf = compute_norm(err,x)
print(N,l1,l2,linf)

		
	
	
	
	
		


