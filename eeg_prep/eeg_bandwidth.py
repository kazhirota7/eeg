from scipy import signal
import scipy.io
import pandas as pd
import numpy as np
import statistics as stat
import matplotlib.pyplot as plt

mat = scipy.io.loadmat('example.mat')

eeg_df = pd.DataFrame(np.vstack(mat['example']))
eeg = list(np.hstack(eeg_df.loc[0,:]))

offset = stat.mean(eeg)
eeg = eeg - offset

plt.plot(eeg)
plt.title("Offset EEG")
plt.show()

fs = 1000
fn = fs / 2


t = []
for i in range(0,len(eeg)):
    t.append(i * (1/fs))
print(np.shape(eeg))
Y = scipy.fft.fft(eeg)
L = np.size(Y)
P2 = abs(Y/L)
P1 = P2[0:L//2+1]
P1[1:] = 2 * P1[1:]
df = fs/len(eeg)
f = []
for i in range(L//2+1):
    f.append(i*df)
plt.plot(f,P1)
plt.xlim(0,100)
plt.title("f vs P1")
plt.show()

# get EEG bandwidth for delta, theta, alpha, beta

filter_order = 3

# delta (1-4Hz)
b, a = signal.butter(filter_order, [1/fn, 4/fn], btype='bandpass')
delta = signal.filtfilt(b, a, eeg)

plt.plot(t, delta)
plt.title("Delta (1-4Hz)")
plt.show()

# theta (4-8Hz)
b, a = signal.butter(filter_order, [4/fn, 8/fn], btype='bandpass')
theta = signal.filtfilt(b, a, eeg)

plt.plot(t, theta)
plt.title("Theta (4-8Hz)")
plt.show()

# alpha (8-12Hz)
b, a = signal.butter(filter_order, [8/fn, 12/fn], btype='bandpass')
alpha = signal.filtfilt(b, a, eeg)

plt.plot(t, alpha)
plt.title("Alpha (8-12Hz)")
plt.show()

# beta (12-35Hz)
b, a = signal.butter(filter_order, [12/fn, 35/fn], btype='bandpass')
beta = signal.filtfilt(b, a, eeg)

plt.plot(t, beta)
plt.title("Beta (12-35Hz)")
plt.show()



