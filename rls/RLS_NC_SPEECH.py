import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf


M = 10
w = np.zeros((M,1))
lamda = 0.9999
delta =4
P = np.identity(M)/delta
#print P

d,fs = sf.read('signal_noise.wav')
v,fs = sf.read('noise.wav')
#print(np.shape(d))
#print(np.shape(v))
#print(fs)

if(len(d) <= len(v)):
	N = len(d)
else:
	N = len(v)
#print(N)

y = np.zeros(N)
e = np.zeros(N)
K = np.zeros((M,1))
#u = np.zeros((filtlen,1))
u = np.concatenate((np.zeros(M-1),v))
print np.shape(u)
for i in range(N):
	x = np.matrix(np.reshape((np.flipud(u[i:i+M])),(M,1)))
	#print np.shape(x)	
	
	K = np.matmul(P,x)/(lamda + np.matmul(np.matmul(x.T,P),x))
	#print K
	
	y[i] = np.matmul(w.T,x)
	#print e[i]
	
	y[i] = d[i] - y[i]
	#print y[i]
	
	w += y[i]*K
	
	P =(P - (np.matmul(np.matmul(K,x.T),P)))/lamda

s_hat = e	
sf.write('output_signal_rls.wav',s_hat,fs)	 
		
	
			 



