import matplotlib.pyplot as plt
import numpy as np
import math

#data_name = 'W20180904131846'
data_name = 'W20180904132846'
#data_name = 'W20180904135846'

data = np.genfromtxt('wind_data/' + data_name + '.csv', delimiter = ',', skip_header = 6001, skip_footer = 6000, usecols = 17)
psdata = np.genfromtxt('gust_data/' + data_name + '_10min.f71', delimiter = [15, 15], skip_header = 1, skip_footer = 1, usecols = 1)
psfdata = np.genfromtxt('gust_data/' + data_name + '_10min.f71', delimiter = [15, 15], skip_header = 1, skip_footer = 1, usecols = 0)

print(psfdata.shape[-1])

sigma_ps = np.sqrt(np.sum(psdata) * 2 * psfdata[1])
print('sigma from spectrum:' + str(sigma_ps))

av_t = np.average(data)
max_t = np.max(data)
std_t = np.std(data)
print('Values from time series')
print(' average: ' + str(av_t))
print(' max: ' + str(max_t))
print(' sigma: ' + str(std_t))
print(' max from average / sigma: ' + str((max_t - av_t) / std_t))

dt = 0.05

time = np.arange(data.shape[-1]) * dt

sp_np = np.fft.rfft(data)
psp_np = np.abs(sp_np) ** 2.0
sigma = np.sqrt(np.sum(psp_np[1:]) * 2.0) / data.shape[-1]
print('sigma:' + str(sigma))
freq = np.fft.rfftfreq(data.shape[-1], d = dt)

print(sigma_ps / sigma)

plt.plot(freq[1:], 2.0 * psp_np[1:] * dt / data.shape[-1], psfdata[1:], psdata[1:])
'''
Relation between numpy power spectrum (Sp_np) v.s. frequency based power spectrum (Sp)
Sp = dt / N * Sp_np * 2 (f > 0)
where dt: sampling time interval, N: number of time series data
'''
plt.yscale('log')
plt.xscale('log')
plt.show()
