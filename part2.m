clc;   
close all;   
clear all;


bit_leng = 500000;
L=2;


SNR_dB = 0:2:20;                     %signal to noise ratio in dB
SNR_dec = 10.^(SNR_dB/10);          %signal to noise ratio in decimal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bit Source%%%%%%%%%%%%%%

bits = randi([0,1],1,bit_leng); %Generate the random bits with length=bit_length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%QPSK Modulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Write out the complex signal for each point on the QPSK constellation. These are the possible signals to transmit
sym_00=exp(1j*pi/4); 
sym_01=exp(1j*3*pi/4); 
sym_11=exp(1j*5*pi/4); 
sym_10=exp(1j*7*pi/4); 

%For any 2 consecutive bits in x, map the bits onto the required qpsk symbol
tx_qpsk_sym = zeros(1,bit_leng/2); %Initialization

j=1;
for i=1:2:length(bits)  
    if bits(i)==0 && bits(i+1)==0
        y=sym_00;
    elseif bits(i)==0 && bits(i+1)==1
        y=sym_01;
    elseif bits(i)==1 && bits(i+1)==1
        y=sym_11;
    elseif bits(i)==1 && bits(i+1)==0
        y=sym_10;
    end
    tx_qpsk_sym(1,j)=y;
    j=j+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rayleigh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 
for i1=1:length(tx_qpsk_sym)
    repeat(1+(i1-1)*L:i1*L)=tx_qpsk_sym(i1);

end
%% fm



ber_sim = zeros(1,length(SNR_dB)); %Initialize the variable to store the BER

for snrloop=1:1:length(SNR_dB)

    noise_signal = (1/sqrt(2))*(randn(1,length(tx_qpsk_sym)*L) + 1j*randn(1,L*length(tx_qpsk_sym))); 
    noise_power = 1/sqrt(SNR_dec(1,snrloop));      %Noise power 
    g = 1/sqrt(2) * (randn(1, L*length(tx_qpsk_sym))) + 1j*randn(1, L*length(tx_qpsk_sym)); 

    rx_signal= g.* repeat + (1/sqrt(2))*(noise_power*noise_signal); %received signal with normalized noise
    for i2=1:length(tx_qpsk_sym)
        demod(i2)=sum(conj(g(1+(i2-1)*L:i2*L)).*rx_signal(1+(i2-1)*L:i2*L))./sum(abs(g(1+(i2-1)*L:i2*L)));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rx_signal=demod;
    %%%%%%%%%%%%%%%%%Maximum likelihood detection%%%%%%%%%%%%%%%%%%%%%%%
    error=[abs(rx_signal-sym_00).',abs(rx_signal-sym_01).',abs(rx_signal-sym_11).',abs(rx_signal-sym_10).'];
    [M,I]=min(error,[],2);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%% Gray demodulation %%%%%%%%%%%%%%%%%%%%%%
    assign_sym=zeros(1,length(rx_signal));
    detec_bits=zeros(1,length(bits));
    j=1;
    for i=1:1:length(rx_signal)
        if I(i,:)==1 %sym_00 is received
            assign_sym(1,i)=sym_00;
            detec_bits(1,j)=0;
            detec_bits(1,j+1)=0;
        elseif I(i,:)==2 %sym_01 is received
            assign_sym(1,i)=sym_01;
            detec_bits(1,j)=0;
            detec_bits(1,j+1)=1;
        elseif I(i,:)==3 %sym_11 is received
            assign_sym(1,i)=sym_11;
            detec_bits(1,j)=1;
            detec_bits(1,j+1)=1;
        elseif I(i,:)==4 %sym_10 is received
            assign_sym(1,i)=sym_10;
            detec_bits(1,j)=1;
            detec_bits(1,j+1)=0;
        end
        j=j+2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ber_sim(1,snrloop) = (sum(bits~=detec_bits))/bit_leng;
%     ber_sim(1,snrloop)=sum(abs(bits-detec_bits))/bit_leng;  
end
    ber_theo=1./(factorial(L)*SNR_dec.^(L)); %Theoretical BER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------Plots---------%
figure;
semilogy(SNR_dB,ber_sim,'^b-',SNR_dB,ber_theo,'dr-'); %Log-based plot
hold on;
grid on;
legend('Simulation','Theory');
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('BER performance of Gray encoding QPSK')







