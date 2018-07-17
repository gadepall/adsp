import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

def lms(N,M,mu,sigma,bh,b):
	J = np.zeros(N)
	J_avg = np.zeros(N)
	e = np.zeros(N)
	y = np.zeros(N)   
	noise = sigma*np.random.randn(len(bh))
	d = bh+noise
	w = np.zeros(M)
	for j in range(M,N):
		u = np.flipud(d[j-M:j])
		y[j] = np.dot(w,u)
		e[j] = b[j-7] - y[j]
		w += mu*e[j] * u
		J[j] += (e[j])**2.0
		
			
	return y,J,w		


N=int(1e6) 
M = 7 # No. of taps
mu = 0.025 # Step size
itr = 50
F = 3.0
snrlen=11
Eb_N0_dB=np.arange(0,snrlen)
sigma=np.sqrt(0.5*(10**(-Eb_N0_dB/10.0)))
h = np.zeros(4)
for n in range(1,len(h)):
	h[n] = 0.5*(1+np.cos(2*np.pi*(n-2)/F))
print h
b=np.zeros(N)
b_hat=np.zeros(N)
simBersoft=np.zeros(snrlen)
nErr=np.zeros(snrlen)

for i in range(snrlen):
	ip=np.random.randint(2,size=N)
	b=2*ip-1
	bh = np.convolve(h,b)
	y,J,w = lms(N,M,mu,sigma[i],bh,b)
	#print w	
	b_hat = 2*(y>=0)-1
	nErr[i] = 0
	for j in range(M,N):
		if b[j-M] != b_hat[j] :
			nErr[i] +=1	


theoryBer=0.5*erfc(np.sqrt((10**(Eb_N0_dB/10.0))))

simBersoft=nErr/float(N)  
#print simBersoft
plt.plot(Eb_N0_dB,theoryBer,'r',Eb_N0_dB,simBersoft,'b')
plt.legend(['theory','Practical'],loc=1)
#plt.plot(Eb_N0_dB,simBersoft,'b')
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for BPSK in raised cosine pulse channel with LMS Adaptive Equalizer ')
plt.grid()
plt.savefig('SNR_Vs_BER.eps')
plt.show()



    
