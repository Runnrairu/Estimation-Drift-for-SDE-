import math
import numpy

def han2(n):
	j=0
	while True:
		n=n/2
		if n==1:
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
	
m=2
k=1
n=2**m+k
l=1000
T=10.0
drift=0.5
sigma2=0.2
X=[0]*(l+1)
Y=[0]*(l+1)
W=[0]*(n)
Z=[0]*(n)
b=[0]*(n)
h= [[0] * (l+1) for i in range(n)]
delta_t=T/l
vol_t=delta_t**(1/2)
for i in range(n):
	m=han2(i)
	k=i-2**m
	for j in range(l+1):
		if j>=k*T/2**m and j<(2*k+1)*T/2**(m+1):
			h[i][j]=(j-k*T/2**m)*((2**m/T)**(1/2))
		elif (2*k+1)*T/2**(m+1)<=j and (k+1)*T/2**m:
			h[i][j]=(T/2**m+2)**(1/2)-((2**m/T)**(1/2))*(j-(2*k+1)+T/2**(m+1))
for i in range (l+1):
	X[i+1]=X[i]+drift*delta_t+sigma*numpy.random.normal(0,vol_t)
for i in range(n):
	W[i]=numpy.random.normal(0,1)
for i in range(n):
	Z[i]=sigma*W[i]/lamda(i)
for i in range(n):
	b[i]=uknaiseki(i)/lamda(i)
fanction()


