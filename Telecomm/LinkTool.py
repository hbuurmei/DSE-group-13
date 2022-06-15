import matplotlib.pyplot as plt
import numpy as np
from scipy import special as sp

## general constants
R_E = 6378  # km
c = 299792458  # m/s
kdB = 1.380649e-23 # m2 kg s-2 K-1

## SC parameters
alt = 1000  # km
f = 3 * 10**9  # Hz, EM wave frequency
L_pr = 0  # [dB] pointing accuracy

## GS parameters
GS_FM = 16  # [dB] Ground station figure of merit

## Parameters of both
L_ldB, L_rdB = 0.7, 0.7 # efficiencies of the transmitter and receiver

## Telecommunication
L_a = -4e-2  # [dB] Worst case atmospheric attenuation
B = 6 * 10**6  # [Hz] signal bandwidth
R = 5 * 10**6  # bit/s

Modulation = np.array([10.6, 11.2, 12.7, 14, 14.5, 18.3, 18.8, 23.3])+5

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
    d = np.sqrt((R+h)**2 - R**2)
    lamb = c/f
    L_s = (lamb/(4 * np.pi * d))**2
    return L_s

L_ldB = to_dB(L_ldB)
L_sdB = to_dB(L_sF(f, alt))
L_rdB = to_dB(L_rdB)
Bdb = to_dB(B)
kdB = to_dB(kdB)
Gt = [6.5]; P_sat = np.arange(0.1, 10, 0.1)
for G_t in Gt:
    PdB = to_dB(P_sat)
    GtdB = G_t
    EIRP = PdB+GtdB+L_ldB
    #SNR_db = EIRP + L_a + L_sdB + L_pr + L_rdB + GS_FM - Bdb - kdB
    EbNo_db = EIRP + L_a + L_sdB + L_pr + L_rdB + GS_FM - to_dB(R) - kdB
    plt.plot(to_normal(PdB), EbNo_db)
plt.grid(True)
plt.xlabel("Power [W]")
plt.ylabel(r"$E_b / N_0$ [dB]")
for m in Modulation:
    plt.plot([P_sat[0], P_sat[-1]], [m, m], linestyle='dashed')
plt.legend(('SSA01', 'BSPK / QPSK / 4-QAM', 'D-BPSK', 'D-QPSK', '8-PSK', '16-QAM', '16-PSK', '64-QAM', '32-PSK'))
plt.show()
#print(SNR_db)
#print(EbNo_db)


# print(np.log10(sp.erfc(np.sqrt(np.arange(1, 15, 1)))))