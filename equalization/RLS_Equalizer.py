import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

def lms(u, d, M, step):

    N = len(u)
    
    # Initialization
    y = np.zeros(N) # Filter output
    e = np.zeros(N) # Error signal
    J = np.zeros(N)
    w = np.zeros((M,1))  # Initialise equaliser
    lamda = 0.9999
    delta = 4
    P = np.identity(M)/delta
    #prnt P
    
    # Equalise
    for i in range(M,N):
		x = np.matrix(np.reshape((np.flipud(u[i-M:i])),(M,1)))
		#print np.shape(x)	
		
		K = np.matmul(P,x)/(lamda + np.matmul(np.matmul(x.T,P),x))
		#print K
		
		y[i] = np.matmul(w.T,x)
		#print e[i]
		
		e[i] = d[i-7] - y[i]
		#print y[i]
		J[i] = (e[i])**2.0
		
		w += e[i]*K   
		
		P =(P - (np.matmul(np.matmul(K,x.T),P)))/lamda

    return y, e, w
    

# Main Program starts from here
N=int(1e6) 
M = 7 # No. of taps
step = 0.003 # Step size

snrlen=11
Eb_N0_dB=np.arange(0,snrlen)
b=np.zeros(N)
b_hat=np.zeros(N)
simBersoft=np.zeros(snrlen)
nErr=np.zeros(snrlen)
F = 3.0
h = np.zeros(3)
for n in range(1,len(h)):
	h[n] = 0.5*(1+np.cos(2*np.pi*(n-1)/F))
#print h

for i in range(snrlen):
    ip=np.random.randint(2,size=N)
    b=2*ip-1
    bh = np.convolve(h,b)

    sigma=np.sqrt(0.5*(10**(-Eb_N0_dB[i]/10.0)))
    #print sigma**2
    
    noise = sigma * np.random.randn(len(bh))
    d = bh+noise
    y, e, w = lms(d,b, M, step)
    b_hat = 2*(y>=0)-1
    nErr[i] = 0
    for j in range(M,N):
		if b[j-M] != b_hat[j] :
			nErr[i] +=1	
	#print nErr[i]		

simBersoft=nErr/float(N)
theoryBer=0.5*erfc(np.sqrt((10**(Eb_N0_dB/10.0))))

print simBersoft
plt.plot(Eb_N0_dB,theoryBer,'r',Eb_N0_dB,simBersoft,'b')
plt.legend(['theory','Practical'],loc=1)
#plt.plot(Eb_N0_dB,simBersoft,'b')
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for BPSK in raised cosine pulse channel with RLS Adaptive Equalizer ')
plt.grid()
plt.savefig('SNR_Vs_BER.eps')
plt.show()



