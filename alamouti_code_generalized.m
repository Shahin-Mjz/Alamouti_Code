

clear; close all; clc

%% General Variables: #####################
M = 2;
n_rx = 2;
coefficient = 1;
iteration = 1e4;
snr_min = 0;
snr_max = 50;
snr = snr_min: snr_max;
%% ########################################

out1 = alamouti_code(M, iteration, coefficient, n_rx);
out2 = alamouti_code(M, iteration , .5^.5, n_rx);

%% Semilogy Plot: ###############################################
semilogy(snr, out1, '-ob', 'DisplayName','Non equal power with MRRC')
hold on 
semilogy(snr, out2, '-or', 'DisplayName', 'Equal power with MRRC')
xlim([0, max(snr)])
ylim([1e-6, 1])
xlabel('SNR (dB)')
ylabel('bit error rate (BER)')
title(['BER for Alamouti Code (1 TX, ', num2str(n_rx), ' RX)'])
grid
legend('FontWeight','bold', 'FontAngle','italic')
%% ##############################################################

function bit_error_rate = alamouti_code(M1, iter, coef, n_rx1)
    snr_min = 0;
    snr_max = 50;
    snr = snr_min: snr_max;
    size_of_symbols = 1e4;   
    ber = zeros(1, length(snr)); ber_sum = zeros(1, length(snr));
    
    for ii = 1: iter
        x = randi([0, M1-1], [size_of_symbols, 1]);   
        s_base = pskmod(x, M1);
        %% Changing Power of Data: #############################
        s_vector = coef * s_base;
        %% #####################################################
        h_vector = (1 / sqrt(2)) * (randn((n_rx1) * size_of_symbols, 1) + 1j*randn((n_rx1) * size_of_symbols, 1)); 
        for snr_index = 1: length(snr)
            num1 = 0;
            s_tild_vector = zeros(size_of_symbols, 1);     
            for i = 0: (size_of_symbols / 2) - 1
                s0 = s_vector(2*i + 1); s1 = s_vector(2*i + 2); 
                s = [s0; s1];
                %%  Creating HA: #############################
                HA = zeros(2*n_rx1, 2);
                num = 0;
                for o = 1: n_rx1
                    num1 = num1 + 2;
                    HA(o + num, 1) = h_vector(num1 - 1); 
                    HA(o + num, 2) = h_vector(num1);
                    HA(o + num + 1, 1) = conj(h_vector(num1));
                    HA(o + num + 1, 2) = -conj(h_vector(num1 - 1));  
                    num = num + 1; 
                end
                %% ###########################################   
                %% Creation of Z: ############################
                z = HA * s; 
                z_vector = awgn(z, snr(snr_index));
                %% ###########################################
                HA_her = HA';
                %% Creating Combining Signals: ###############
                s_tild = HA_her * z_vector;
                s_tild_vector((2*i + 1): (2*i + 2)) = s_tild;
                %% ###########################################
            end
            %% s_hat_vector Is Output of Demodulator:######
            s_hat_vector = pskdemod(s_tild_vector, M1);
            %% ############################################
            %% Bit Error Rate Calculaiton:
            input_bit_data = de2bi(x);
            output_bit_data = de2bi(s_hat_vector);
            number_of_error = 0;
            for f = 1: size_of_symbols
                a = find(output_bit_data(f, :) ~= input_bit_data(f, :));
                number_of_error = number_of_error + length(a);
            end
            ber(snr_index) = number_of_error / (size_of_symbols * log2(M1));
        end
        ber_sum = ber_sum + ber;
    end
    bit_error_rate = ber_sum / iter;
    %% #####################################################
end













