import numpy as np

## general constants
R_E = 6378  # km
c = 299792458  # m/s
k = 1.380649e-23 # m2 kg s-2 K-1
## SC parameters
alt = 1000  # km
f = 3 * 10**9  # Hz, EM wave frequency
P_sat = 1  # Power available to the TT&C subsystem
G_t = 2
L_pr = 0  # [dB] pointing accuracy

## GS parameters
GS_FM = 16  # [dB] Ground station figure of merit

## Parameters of both
L_l, L_r = 0.7, 0.7 # efficiencies of the transmitter and receiver

## Telecommunication
L_a = -4e-2  # [dB] Worst case atmospheric attenuation
B = 6 * 10**6  # [Hz] signal bandwidth
D = 5 * 10**5  # bits
T_comm = 1 # s
R = D/T_comm # bit/s

to_dB = lambda x: 10 * np.log10(x)
to_normal = lambda x: 10**(x/10)

def L_sF(f, h=1000, R=R_E):
    '''
    @:param f, the frequency in Hz
    @:param h, the altitude in km
    @:param R, the radius of the planet
    '''
    h *= 1000
    R *= 1000
    d = np.sqrt((R+h)**2 -R**2)
    lamb = c/f
    L_s = (lamb/(4 * np.pi * d))**2
    return L_s

P = to_dB(P_sat)
L_l = to_dB(L_l)
L_s = to_dB(L_sF(f, alt))
L_r = to_dB(L_r)
B = to_dB(B)
k = to_dB(k)


Bguess, err = 6e6, 10
B1 = Bguess
while err > 1e-4:
    SNR = 2**(R*0.5 /B1)
    SNR_db = to_dB(SNR)
    BdB = P + L_l + L_a + L_s + + L_pr + L_r + GS_FM - SNR_db - k
    B2 = 10**(BdB/10)
    err = abs(B2-B1)
    print(B1, SNR)
    B1 = B2

EbNo = SNR + (B - R)
SNR_db = P + L_l + L_a + L_s + + L_pr + L_r + GS_FM - BdB - k
EbNo_db = P + L_l + L_a + L_s + + L_pr + L_r + GS_FM - to_dB(R) - k

print(EbNo_db)