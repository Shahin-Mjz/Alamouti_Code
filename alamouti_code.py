import numpy as np
import komm
from math import sqrt, log2
from matplotlib import pyplot as plt

## General Variables: ######################################
M = 2
n_rx = 2
coefficent = 1
iteration = int(1e3)
snr_min, snr_max = 0, 50
snr_list = range(snr_min, snr_max + 1)
snr = np.array(snr_list).reshape((1, len(snr_list)))
############################################################


def alamouti_code(M1, iter, coef, n_rx1, snr_min1, snr_max1):
    snr1 =  range(snr_min1, snr_max1 + 1)
    size_of_data = int(1e3)
    ber, ber_sum = np.zeros((1, len(snr1)), dtype=float), np.zeros((1, len(snr1)), dtype=float)
    for ii in range(1, iter + 1):
        x = np.random.randint(0, M1, size = (size_of_data, 1))
        psk = komm.PSKModulation(M1)
        s_base = psk.modulate(x)
        ## Changing Power of Data: ############################################################################################################
        s_vector1 = coef * s_base
        #######################################################################################################################################
        s_vector = s_vector1.reshape((len(s_vector1), 1))
        h_vector = 1/sqrt(2) * ((np.random.normal(size=(n_rx1*size_of_data, 1))) + 1j*(np.random.normal(size=(n_rx1*size_of_data, 1))))
        for snr_value in range(snr_min1, snr_max1 + 1):
            num1 = 0
            s_tild_vector = np.zeros((size_of_data, 1), dtype=complex)
            for i in range(0, int(size_of_data / 2)):
                s0, s1 = s_vector[2*i, 0], s_vector[2*i + 1, 0]
                s = np.array([[s0], [s1]]) 
                ## Creating HA: ########################################################
                HA = np.zeros((2*n_rx1, 2), dtype=complex)
                num2 = 0
                for o in range(0, n_rx1):
                    num1 = num1 + 2
                    HA[o + num2, 0] = h_vector[num1 - 2, 0]
                    HA[o + num2, 1] = h_vector[num1 - 1, 0]
                    HA[o + num2 + 1, 0] = np.conj(h_vector[num1 - 1, 0])
                    HA[o + num2 + 1, 1] = -np.conj(h_vector[num1 - 2, 0])
                    num2 += 1
                ########################################################################
                ## Creating z:###########################################
                z = HA @ s
                SNR_linear = 10 ** (snr_value / 10)
                awgn = komm.AWGNChannel(snr=SNR_linear)  
                z_vector = awgn(z)
                #########################################################
                HA_hermitian = HA.T.conj()
                ## Creating Combining Signals: ##########################
                s_tild = HA_hermitian @ z_vector
                s_tild_vector[2*i, 0] = s_tild[0, 0]
                s_tild_vector[2*i + 1, 0] = s_tild[1, 0]
                #########################################################
            ## s_hat_vector Is Output of demodulator: ######################
            s_hat_vector1 = psk.demodulate(s_tild_vector)
            s_hat_vector = s_hat_vector1.reshape((len(s_hat_vector1), 1))
            ################################################################
            ## Bit Error Rate Calculation: ######################################
            number_of_errors = 0
            for p in range(size_of_data):
                for q in range(len(bin(x[p, 0])[2:])):
                    input_bit_data = int(bin(x[p, 0])[2 + q])
                    output_bit_data = int(bin(s_hat_vector[p, 0])[2 + q])
                    if input_bit_data != output_bit_data:
                        number_of_errors += 1
            snr_index = snr_value
            ber[0, snr_index] = number_of_errors / (size_of_data * log2(M1))
        ber_sum += ber
    bit_error_rate = ber_sum / iter
    #############################################################################
    return bit_error_rate
    

out = alamouti_code(M, iteration, coefficent, n_rx, snr_min, snr_max)

## Semilogy Plot: ########################################################
fig, ax = plt.subplots()
ax.semilogy(snr, out, marker = '.', color='blue')
ax.set_title('Alamouti Code for BPSK Modulation and n_rx = ' + str(n_rx))
ax.set_xlabel('SNR(dB)')
ax.set_ylabel('Bit error rate(BER)')
ax.set_xlim([snr_min + 1, snr_max])
ax.set_ylim([1e-6, 1])
plt.grid(True, which='both')
plt.show()
##########################################################################