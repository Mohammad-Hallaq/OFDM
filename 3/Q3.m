clc
close all
clear

%TRANSMITTER

N_bps = 64;
M_ord = 2;
N_scu = 4;
T_N_b = 256;
Block_size = 16;
n = 7;
k = 4;
M = 2;
K = log2(M);
Len = 1e4;
alphabet = 0 : M - 1; 
A = 1; 
teta0 = 0; 
Ts = 1/1e5;
Fd = 500;

CyclicPrefix_size = ceil(Block_size/10);

data = randsrc(T_N_b,1,[0:M_ord-1]);
figure(1);
stem(data);
enc = encode(data, n,k,'hamming/binary'); 
ss = 2*enc - 1;   
int_1=reshape(ss,[],size(ss,1));
int_2=transpose(int_1);
int_3=reshape(int_2,1,[]);        
sym = int_3*sqrt(k/n);


resahped_data = reshape(sym,[7*Block_size, N_scu]);

sub_data_1 = resahped_data(:,1);
sub_data_2 = resahped_data(:,2);
sub_data_3 = resahped_data(:,3);
sub_data_4 = resahped_data(:,4);

figure(3);
subplot(2,2,1); stem(sub_data_1);
subplot(2,2,2); stem(sub_data_2);
subplot(2,2,3); stem(sub_data_3);
subplot(2,2,4); stem(sub_data_4);

spect_sub_data_1 = ifft(sub_data_1);
spect_sub_data_2 = ifft(sub_data_2);
spect_sub_data_3 = ifft(sub_data_3);
spect_sub_data_4 = ifft(sub_data_4);

figure(4);
subplot(2,2,1); plot(real(spect_sub_data_1));
subplot(2,2,2); plot(real(spect_sub_data_2));
subplot(2,2,3); plot(real(spect_sub_data_3));
subplot(2,2,4); plot(real(spect_sub_data_4));

prefix_1 = spect_sub_data_1(length(spect_sub_data_1) - CyclicPrefix_size +1 :length(spect_sub_data_1));
prefix_2 = spect_sub_data_2(length(spect_sub_data_2) - CyclicPrefix_size +1 :length(spect_sub_data_2));
prefix_3 = spect_sub_data_3(length(spect_sub_data_3) - CyclicPrefix_size +1 :length(spect_sub_data_3));
prefix_4 = spect_sub_data_4(length(spect_sub_data_4) - CyclicPrefix_size +1 :length(spect_sub_data_4));

ofdm_block_1 = vertcat(prefix_1, spect_sub_data_1);
ofdm_block_2 = vertcat(prefix_2, spect_sub_data_2);
ofdm_block_3 = vertcat(prefix_3, spect_sub_data_3);
ofdm_block_4 = vertcat(prefix_4, spect_sub_data_4);

figure(5);
subplot(2,2,1); plot(real(ofdm_block_1));
subplot(2,2,2); plot(real(ofdm_block_2));
subplot(2,2,3); plot(real(ofdm_block_3));
subplot(2,2,4); plot(real(ofdm_block_4));

OFDM_signal = vertcat(ofdm_block_1, ofdm_block_2, ofdm_block_3,ofdm_block_4);
figure(6);
stem(OFDM_signal);
%%
close all
%CHANNEL
  
x =zeros(length(OFDM_signal) ,1);
 
SNR_mat = -4:20;
received_signal_mat = zeros(length(OFDM_signal), length(SNR_mat));
 
 for i=1:length(SNR_mat)
     
     noise = awgn(x,SNR_mat(i));
     received_signal_mat(:,i) = OFDM_signal + noise ;
     
 end


%%

%RECEIVER

for j=1:length(SNR_mat)
    
    r = received_signal_mat(:, j);
    
    data_out_1 = r(1: 7*Block_size + CyclicPrefix_size);
    
    data_out_2 = r(7*Block_size + CyclicPrefix_size + 1 : 2*(7*Block_size + CyclicPrefix_size));
    
    data_out_3 = r(2*(7*Block_size + CyclicPrefix_size)+1 : 3*(7*Block_size + CyclicPrefix_size));
    
    data_out_4 = r(3*(7*Block_size + CyclicPrefix_size)+1 : 4*(7*Block_size + CyclicPrefix_size));
    
    data_block_out_1 = data_out_1(CyclicPrefix_size + 1 : 7*Block_size + CyclicPrefix_size);
    data_block_out_2 = data_out_2(CyclicPrefix_size + 1 : 7*Block_size + CyclicPrefix_size);
    data_block_out_3 = data_out_3(CyclicPrefix_size + 1 : 7*Block_size + CyclicPrefix_size);
    data_block_out_4 = data_out_4(CyclicPrefix_size + 1 : 7*Block_size + CyclicPrefix_size);
    
    inv_data_out_1 = fft(data_block_out_1);
    inv_data_out_2 = fft(data_block_out_2);
    inv_data_out_3 = fft(data_block_out_3);
    inv_data_out_4 = fft(data_block_out_4);
    
    data_serial_out = vertcat(inv_data_out_1, inv_data_out_2, inv_data_out_3, inv_data_out_4);
    
     recieved_data1 = (sign(data_serial_out)+1)/2;
     
     recieved_data = pskdemod(recieved_data1, M_ord);
     
     de_int_1=reshape(recieved_data,size(recieved_data,1),[]);
     de_int_2=transpose(de_int_1);
     de_int_3=reshape(de_int_2,[],1);
     
     dec = decode(de_int_3, n,k, 'hamming/binary');
     

    [num(j), ratio(j)] = biterr(data, dec);
    
end



plot(SNR_mat, ratio);

        
        

