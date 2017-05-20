


import math
import numpy

def han2(n):
	j=0
	while True:
		n=n/2
		if n==1 or n==0:
			break
		else:
			j=j+1
	return j

def lamda(i):
	return (sigma**2)*((1/T)+T**2)/(3*(2**(2*han2(i)+3)))
			

def uknaiseki(i):
	integralL2=0
	for j in range(l):
		integralL2=integralL2+h[i][j+1]*X[j+1]*delta_t
	return integralL2

m=6
k=20
n=2**m+k
l=1000
T=10.0
drift=0.5
sigma=0.2
X=[0]*(l+2)
Y=[0]*(l+2)
W=[0]*(n)
Z=[0]*(n)
b=[0]*(n)
h= [[0] * (l+2) for i in range(n)]
lam=[0]*n
delta_t=T/l
vol_t=delta_t**(1/2)

for i in range(n):
	m=han2(i)
	k=i-2**m
	for j in range(l+1):
		t=float(j/l)
		if t>=k*T/2**m and t<(2*k+1)*T/2**(m+1):
			h[i][j]=(t-k*T/2.0**m)*((2**m/T)**(0.5))
		elif (2*k+1)*T/2**(m+1)<=t and (k+1)*T/2**m>t:
			h[i][j]=(T/2**m+2)**(1/2.0)-((2**m/T)**(0.5))*(t-(2*k+1)+T/2**(m+1))
for i in range (l+1):
	X[i+1]=X[i]+drift*delta_t+sigma*numpy.random.normal(0,vol_t)
for i in range(n):
	W=numpy.random.normal(0,1,n)
for i in range(n):
	lam[i]=lamda(i)
for i in range(n):
	Z[i]=sigma*W[i]/lam[i]
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
		de_sum=de_sum+gx**((2-n)/2.0-1)*(b[i]+2*lam[i])*sigma**2*h[i][j]/lam[i]
	malliavin_delive[j]=de_sum/gx**((2-n)/2)
	drift_estimator[j]=malliavin_delive[j]+X[j]
print drift_estimator

