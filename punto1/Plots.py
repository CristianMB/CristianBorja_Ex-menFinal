import numpy as np
import matplotlib.pyplot as plt

MCMC_1 = np.loadtxt('sample_0.dat')
MCMC_2 = np.loadtxt('sample_1.dat')
MCMC_3 = np.loadtxt('sample_2.dat')
MCMC_4 = np.loadtxt('sample_3.dat')
MCMC_5 = np.loadtxt('sample_4.dat')
MCMC_6 = np.loadtxt('sample_5.dat')
MCMC_7 = np.loadtxt('sample_6.dat')
MCMC_8 = np.loadtxt('sample_7.dat')

MCMCs = [MCMC_1, MCMC_2, MCMC_3, MCMC_4, MCMC_5, MCMC_6, MCMC_7, MCMC_8]  

mu = 0.0
sigma = 1.0
x = np.arange(-5, 5, 0.1)

plt.figure()
plt.title('Sampling')
plt.plot(x, np.exp(-(x-mu)*2/(2*sigma*2))/(sigma*np.sqrt(2*np.pi)), linewidth=2, color='darkgreen')
for i in MCMCs:
    plt.hist(i, bins=100, normed=True, color='blue')
plt.savefig("MCMC_Sampling.pdf")
