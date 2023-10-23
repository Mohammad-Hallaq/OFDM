function [ber, numBits] = exp1_MonteCarlo2(EbNo, maxNumErrs, maxNumBits)
% Import Java class for BERTool.
import com.mathworks.toolbox.comm.BERTool;

% Initialize variables related to exit criteria.
berVec = zeros(3,1); % Updated BER values

N_bps = 64;
M_ord = 4;
N_scu = 8;
T_N_b = 8192;
Block_size = 32;
CyclicPrefix_size = ceil(Block_size/10);


%Len = 64;

%g = ones(1,fs)/sqrt(fs);
% --- Set up parameters. ---
% --- INSERT YOUR CODE HERE.
% Simulate until number of errors exceeds maxNumErrs
% or number of bits processed exceeds maxNumBits.
while((berVec(2) < maxNumErrs) && (berVec(3) < maxNumBits))
        data = randsrc(T_N_b,1,[0:M_ord-1]); % information bits
        QPSK_data = pskmod(data, M_ord); % QPSK modulation
        resahped_data = reshape(QPSK_data,[32*Block_size, N_scu]);
        sub_data_1 = resahped_data(:,1);
        sub_data_2 = resahped_data(:,2);
        sub_data_3 = resahped_data(:,3);
        sub_data_4 = resahped_data(:,4);
        sub_data_5 = resahped_data(:,5);
        sub_data_6 = resahped_data(:,6);
        sub_data_7 = resahped_data(:,7);
        sub_data_8 = resahped_data(:,8);
        
        spect_sub_data_1 = ifft(sub_data_1);
        spect_sub_data_2 = ifft(sub_data_2);
        spect_sub_data_3 = ifft(sub_data_3);
        spect_sub_data_4 = ifft(sub_data_4);
        spect_sub_data_5 = ifft(sub_data_5);
        spect_sub_data_6 = ifft(sub_data_6);
        spect_sub_data_7 = ifft(sub_data_7);
        spect_sub_data_8 = ifft(sub_data_8);
        
        prefix_1 = spect_sub_data_1(length(spect_sub_data_1) - CyclicPrefix_size +1 :length(spect_sub_data_1));
        prefix_2 = spect_sub_data_2(length(spect_sub_data_2) - CyclicPrefix_size +1 :length(spect_sub_data_2));
        prefix_3 = spect_sub_data_3(length(spect_sub_data_3) - CyclicPrefix_size +1 :length(spect_sub_data_3));
        prefix_4 = spect_sub_data_4(length(spect_sub_data_4) - CyclicPrefix_size +1 :length(spect_sub_data_4));
        prefix_5 = spect_sub_data_5(length(spect_sub_data_5) - CyclicPrefix_size +1 :length(spect_sub_data_5));
        prefix_6 = spect_sub_data_6(length(spect_sub_data_6) - CyclicPrefix_size +1 :length(spect_sub_data_6));
        prefix_7 = spect_sub_data_7(length(spect_sub_data_7) - CyclicPrefix_size +1 :length(spect_sub_data_7));
        prefix_8 = spect_sub_data_8(length(spect_sub_data_8) - CyclicPrefix_size +1 :length(spect_sub_data_8));

        ofdm_block_1 = vertcat(prefix_1, spect_sub_data_1);
        ofdm_block_2 = vertcat(prefix_2, spect_sub_data_2);
        ofdm_block_3 = vertcat(prefix_3, spect_sub_data_3);
        ofdm_block_4 = vertcat(prefix_4, spect_sub_data_4);
        ofdm_block_5 = vertcat(prefix_5, spect_sub_data_5);
        ofdm_block_6 = vertcat(prefix_6, spect_sub_data_6);
        ofdm_block_7 = vertcat(prefix_7, spect_sub_data_7);
        ofdm_block_8 = vertcat(prefix_8, spect_sub_data_8);
        
        OFDM_signal = vertcat(ofdm_block_1, ofdm_block_2, ofdm_block_3,ofdm_block_4,ofdm_block_5,ofdm_block_6,ofdm_block_7,ofdm_block_8);
        
        %channel = randn([2,1])+ 1i*randn([2,1]);
        %filtered_ofdm = filter(channel,1,OFDM_signal);
        x =zeros(length(OFDM_signal) ,1);
        noise = awgn(x, EbNo + 3);
        received_signal_mat = OFDM_signal + noise ;
        
        
        r = received_signal_mat;
        data_out_1 = r(1: 32*Block_size + CyclicPrefix_size);
        data_out_2 = r(32*Block_size + CyclicPrefix_size + 1 : 2*(32*Block_size + CyclicPrefix_size));
        data_out_3 = r(2*(32*Block_size + CyclicPrefix_size)+1 : 3*(32*Block_size + CyclicPrefix_size));
        data_out_4 = r(3*(32*Block_size + CyclicPrefix_size)+1 : 4*(32*Block_size + CyclicPrefix_size));
        data_out_5 = r(4*(32*Block_size + CyclicPrefix_size)+1 : 5*(32*Block_size + CyclicPrefix_size));
        data_out_6 = r(5*(32*Block_size + CyclicPrefix_size)+1 : 6*(32*Block_size + CyclicPrefix_size));
        data_out_7 = r(6*(32*Block_size + CyclicPrefix_size)+1 : 7*(32*Block_size + CyclicPrefix_size));
        data_out_8 = r(7*(32*Block_size + CyclicPrefix_size)+1 : 8*(32*Block_size + CyclicPrefix_size));

        data_block_out_1 = data_out_1(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_2 = data_out_2(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_3 = data_out_3(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_4 = data_out_4(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_5 = data_out_5(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_6 = data_out_6(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_7 = data_out_7(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);
        data_block_out_8 = data_out_8(CyclicPrefix_size + 1 : 32*Block_size + CyclicPrefix_size);

        inv_data_out_1 = fft(data_block_out_1);
        inv_data_out_2 = fft(data_block_out_2);
        inv_data_out_3 = fft(data_block_out_3);
        inv_data_out_4 = fft(data_block_out_4);
        inv_data_out_5 = fft(data_block_out_5);
        inv_data_out_6 = fft(data_block_out_6);
        inv_data_out_7 = fft(data_block_out_7);
        inv_data_out_8 = fft(data_block_out_8);

        data_serial_out = vertcat(inv_data_out_1, inv_data_out_2, inv_data_out_3, inv_data_out_4, inv_data_out_5, inv_data_out_6, inv_data_out_7, inv_data_out_8);
        demod_serial_out = pskdemod(data_serial_out, M_ord);
        [num, ratio] = biterr(data, demod_serial_out);
        
      
        berVec(2) = berVec(2) + num;
        berVec(3) = berVec(3) + T_N_b;
   % Check if the user clicked the Stop button of BERTool.
   if (BERTool.getSimulationStop)
      break;
   end

   % --- Proceed with simulation.
   % --- Be sure to update totErr and numBits.
   % --- INSERT YOUR CODE HERE.
end % End of loop

berVec(1) = berVec(2)/berVec(3);

% Assign values to the output variables.
ber = ratio;
numBits = berVec(3);