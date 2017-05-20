import math
import numpy

def han2(n):
	j=1
	while True:
		
		if 2**j>n:
			break
		else:
			j=j+1
	return j-1

def lamda(i):
	return sigma*T/(math.pi*(n-0.5))
			

def uknaiseki(i):
	integralL2=0
	for j in range(l):
		integralL2=integralL2+h[i][j+1]*X[j+1]*delta_t
	return integralL2


n=100
l=1000
T=10.0
drift=0.5
sigma=1
X=[0]*(l+2)
W=[0]*(n)
Z=[0]*(n)
b=[0]*(n)
h= [[0] * (l+2) for i in range(n+1)]
lam=[0]*n
delta_t=T/l
vol_t=delta_t**0.5

for i in range(n):
	m=han2(i)
	for j in range(l+1):
		t=j*T/l
		h[i][j]=(2*T)**0.5*math.sin((m-0.5)*math.pi*t/T)/(sigma*math.pi*(m-0.5))
			

for i in range (l+1):
	X[i+1]=X[i]+drift*delta_t+sigma*numpy.random.normal(0,vol_t)
for i in range(l+1):
	W2[i]=numpy.random.normal(0,1)
for i in range(n):
	W[i]=numpy.random.normal(0,1)
for i in range(n):
	lam[i]=lamda(i)
for i in range(n):
	Z[i]=sigma**2*W[i]/lam[i]
for i in range(n):
	b[i]=uknaiseki(i)/lam[i]

gx=0

for i in range(n):
	gx=gx+(b[i]+Z[i])**2

malliavin_delive=[0]*(l+1)
drift_estimator=[0]*(l+1)
for j in range(l+1):
	de_sum=0
	for i in range(n):
		de_sum=de_sum+((2-n)/2)*gx**((2-n)/2.0-1)*(b[i]+2*Z[i])*(sigma**2)*h[i][j]/lam[i]
	malliavin_delive[j]=de_sum/gx**((2-n)/2)
	drift_estimator[j]=malliavin_delive[j]+X[j]
print X
print malliavin_delive
print drift_estimator
