import numpy as np

i = 0
M = 3
Z,W = np.polynomial.legendre.leggauss(M)

f = lambda z,n : 3*z+np.exp(-n)

x = lambda z,n : 0.5 - z/2
y = lambda z,n : 1/4*(1-z)*(1+n)

for wz,z in zip(W,Z):
	for wn,n in zip(W,Z):
		i+=((1-z)/8)*f(x(z,n),y(z,n))*wz*wn
print(i,1/2+np.exp(-1))
# n = 3
# z,w = np.polynomial.legendre.leggauss(n)
# c = []
# x = []
# y = []
# for k in range(n):
# 	ck = 0
# 	xk = 0
# 	yk = 0
# 	for i in range(n):
# 		for j in range(n):
# 			ck += (1-z[i])/8*w[i]*w[j]
# 			xk += (1+z[i])/2
# 			yk += (1-z[i])*(1+z[j])/4
# 	c.append(ck)
# 	x.append(xk)
# 	y.append(yk)

# print(c,'\n',x,'\n',y)