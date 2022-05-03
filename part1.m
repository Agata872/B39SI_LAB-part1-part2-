clc;   
close all;   
clear all;

fm = 91;
b = 0.5;
N1 = 9;
N2 = 10;
fs = 270.8e3;
ts = 1/fs;
bit_leng = 500000;
Rb = 541.6e3;


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
c1 = sqrt((2*b)/N1);   
f1 = zeros(1,N1);       
theta1 = 2*pi*rand(1,N1);
for n1 = 1:1:N1
    f1(1,n1) = fm*sin((pi*(n1-(1/2)))/(2*N1));
end

c2 = sqrt ((2*b)/N2); 
f2 = zeros(1,N2);       
theta2 = 2*pi*rand(1,N2);
for n2 = 1:1:N2
    f2(1,n2)  = fm*sin((pi*(n2-(1/2)))/(2*N2));
end

T_sim = bit_leng / Rb;
t = 0:ts:(T_sim-ts);

gn_1=zeros(size(t,2),size(f1,2));  %size(A,dim) returns the length of dimension dim when dim is a positive integer scalar. 
gn_2=zeros(size(t,2),size(f2,2));
for index=1:1:size(t,2) %The number of rows.
    gn_1(index,:) = c1.*cos((2.*pi.*f1.*t(1,index))+theta1);
    gn_2(index,:) = c2.*cos((2.*pi.*f2.*t(1,index))+theta2);
end
g1=sum(gn_1,2); % sum(A,dim) returns the sum along dimension dim. 
g2=sum(gn_2,2);
g = g1 + 1i.*g2; 
alpha = (abs(g));

ber_sim = zeros(1,length(SNR_dB)); %Initialize the variable to store the BER

for snrloop=1:1:length(SNR_dB)

    noise_signal = (1/sqrt(2))*(randn(1,length(tx_qpsk_sym)) + 1j*randn(1,length(tx_qpsk_sym))); 
    noise_power = 1/sqrt(SNR_dec(1,snrloop));      %Noise power 
    rx_signal= alpha.'.* tx_qpsk_sym + (1/sqrt(2))*(noise_power*noise_signal); %received signal with normalized noise
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    ber_theo=(1/2)*(1-sqrt(2*b*SNR_dec./(1+2*b*SNR_dec))); %Theoretical BER
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







