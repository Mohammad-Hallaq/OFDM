clearvars; close all; clc;
M = 4;
T_N_b = 8192;
N_scu = 4;
k = log2(M);
N_sym = T_N_b/k;
Block_size = 16;
S_P_B = Block_size/k;
blocks_per_sub = N_sym/N_scu/S_P_B;
symbols_per_sub = blocks_per_sub*S_P_B;
cp = ceil(Block_size/10);
alphabet = 0:M-1;
S_PSK = exp(1i * 2 * pi * alphabet / M);
dist = zeros(1,M);
EbN0 = 1:20;

for h = 1:length(EbN0)
    numErrs = 0;
    numBits = 0;
while((numErrs < 50) && (numBits < 1e12))

    data = randsrc(1,N_sym,alphabet);
    data = data';
    modulated_data= pskmod(data,M);
    S2P = reshape(modulated_data,symbols_per_sub,N_scu);

    for i=1:N_scu
        totifft(:,i) = ifft(S2P(:,i));
    end
    y = totifft(:);
    sent_par = reshape(y,S_P_B,blocks_per_sub,N_scu);
    sent_par_cp = cat(1,sent_par,sent_par(1:cp,:,:));
    sent_ser_cp = sent_par_cp(:);

    train_len = 200;
    data_len = length(sent_ser_cp);
    train = randsrc(1,train_len,alphabet);
    train = train';
    sent_frame = [train;sent_ser_cp];
    channel_B = [0.407 0.815 0.407];
    delayB = (length(channel_B)-1)/2;
    eff_signalB = conv(channel_B, sent_frame);
    noised_data = awgn(eff_signalB,EbN0(h) + 10*log10(k));

    recieved_data = noised_data(delayB+1:length(noised_data)-delayB);
    noised_train = recieved_data(1:train_len);
    recieved_without_train = recieved_data(train_len+1:end);
    recieved_par_cp = reshape(recieved_without_train,S_P_B + cp,blocks_per_sub,N_scu);
    recieved_par = recieved_par_cp(1:S_P_B,:,:);
    ioi = reshape(recieved_par,symbols_per_sub,N_scu);
    for i=1:N_scu
        fft_data(:,i) = fft(ioi(:,i));
    end
    recieved_ser = fft_data(:);

    N = 2;
    delta = 0.016;
    recieved_frame = [noised_train;recieved_ser];
    c = zeros(2*N,length(recieved_frame));
    c(:,1) = eps*ones(1,2*N);
    for i=1:length(recieved_frame)-2*N
        yk = sum(recieved_frame(i:i+2*N-1).*c(:,i));
        if(i<train_len)
            ek = sent_frame(i)-yk;
        else
            for j = 1 : M
                dist(j) = abs(yk - S_PSK(j));
            end
            [dis, ind] = min(dist);
            ek = dis;
            demoulated_data = alphabet(ind);
            nErrors = biterr(demoulated_data, data(i-train_len+1), k);
            numErrs = numErrs + nErrors;
            numBits = numBits + T_N_b;
        end
        c(:,i+1) = c(:,i)+delta*ek.*noised_data(i:i+2*N-1);
    end
end 
ber(h) = numErrs / numBits;
end
figure(3);
plot(EbN0,ber);
grid on;