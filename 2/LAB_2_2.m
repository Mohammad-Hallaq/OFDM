clc
close all
clear

%TRANSMITTER

N_bps = 64;
M_ord = 4;
N_scu = 8;
T_N_b = 8192;
Block_size = 32;

CyclicPrefix_size = ceil(Block_size/10);

data = randsrc(T_N_b,1,[0:M_ord-1]);
figure(1);
stem(data);

QPSK_data = pskmod(data, M_ord);
figure(2);
stem(QPSK_data);

resahped_data = reshape(QPSK_data,[32*Block_size, N_scu]);

sub_data_1 = resahped_data(:,1);
sub_data_2 = resahped_data(:,2);
sub_data_3 = resahped_data(:,3);
sub_data_4 = resahped_data(:,4);
sub_data_5 = resahped_data(:,5);
sub_data_6 = resahped_data(:,6);
sub_data_7 = resahped_data(:,7);
sub_data_8 = resahped_data(:,8);


figure(3);
subplot(2,2,1); stem(sub_data_1);
subplot(2,2,2); stem(sub_data_2);
subplot(2,2,3); stem(sub_data_3);
subplot(2,2,4); stem(sub_data_4);

spect_sub_data_1 = ifft(sub_data_1);
spect_sub_data_2 = ifft(sub_data_2);
spect_sub_data_3 = ifft(sub_data_3);
spect_sub_data_4 = ifft(sub_data_4);
spect_sub_data_5 = ifft(sub_data_5);
spect_sub_data_6 = ifft(sub_data_6);
spect_sub_data_7 = ifft(sub_data_7);
spect_sub_data_8 = ifft(sub_data_8);




figure(4);
subplot(2,2,1); plot(real(spect_sub_data_1));
subplot(2,2,2); plot(real(spect_sub_data_2));
subplot(2,2,3); plot(real(spect_sub_data_3));
subplot(2,2,4); plot(real(spect_sub_data_4));

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

figure(5);
subplot(2,2,1); plot(real(ofdm_block_1));
subplot(2,2,2); plot(real(ofdm_block_2));
subplot(2,2,3); plot(real(ofdm_block_3));
subplot(2,2,4); plot(real(ofdm_block_4));

OFDM_signal = vertcat(ofdm_block_1, ofdm_block_2, ofdm_block_3,ofdm_block_4,ofdm_block_5,ofdm_block_6,ofdm_block_7,ofdm_block_8);
figure(6);
stem(OFDM_signal);
%%
close all
%CHANNEL

channel_B=[0.407 0.815 0.407];
filtered_ofdm = conv(channel_B, OFDM_signal);

x =zeros(length(filtered_ofdm) ,1);



SNR_mat = -4:20;
received_signal_mat = zeros(length(filtered_ofdm), length(SNR_mat));
 
 for i=1:length(SNR_mat)
     
     noise = awgn(x,SNR_mat(i));
     received_signal_mat(:,i) = filtered_ofdm + noise ;
     
 end
 
 

%%

clc ;

ebno=10;
N=1;
delta=0.01;
train_Len=200;
data_Len=1e5;
train=randsrc(1,train_Len);
data= OFDM_signal; 
train_data=[train data'];
channel_B=[0.407 0.815 0.407];
delay=(length(channel_B)-1)/2;
eff_signal=conv(channel_B, train_data);
recive_signal=awgn(eff_signal,ebno+3);
X=zeros(1,data_Len);
c=eps*ones(2*N+1,length(train_data));
for i=1:length(train_data)-2*delay
    yk=c(:,i)'*recive_signal(i:2*delay+i)';
    zk=sign(yk);
if(i<train_Len+1)
  ek=train(i)-yk;
else
   ek=zk-yk;
    X(i-train_Len)=yk;
end
c(:,i+1)=c(:,i)+ek*delta*recive_signal(i:2*delay+i)';
end
compare=data-sign(X);
error=nnz(compare);
err_rate=error/length(compare);
hist(X,100);
grid on
title('Recieved data histogram channel B with the Equalizer');
%%
%RECEIVER

X=zeros(1,length(filtered_ofdm));
ebno=10;
N=1;
delta=0.01;
train_Len=200;
data_Len=1e5;
train=randsrc(1,train_Len);
data=randsrc(1,data_Len);
train_data=[train data];
delay=(length(channel_B)-1)/2;
c=eps*ones(2*N+1,length(train_data));
for i=1:length(train_data)-2*delay
    yk=c(:,i)'*received_signal_mat(i:2*delay+i)';
    zk=sign(yk);
if(i<train_Len+1)
  ek=train(i)-yk;
else
   ek=zk-yk;
    X(i-train_Len)=yk;
end
c(:,i+1)=c(:,i)+ek*delta*received_signal_mat(i:2*delay+i)';
end

for j=1:length(SNR_mat)
    
    r = received_signal_mat(:, j);
    
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
    
    [num(j), ratio(j)] = biterr(data, demod_serial_out);
    
end



plot(SNR_mat, ratio);








