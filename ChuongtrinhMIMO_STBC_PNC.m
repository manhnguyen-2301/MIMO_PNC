%===============================================================

clear all
rand('seed',0);
randn('seed',0);

%>>>>> SIMULATION SETUP <<<<<
% Simulation parameters
num_node = 2;           % number of end_nodes
ant_node = 2;           % number of antenna of a node
num_relay = 1;           % number of relay
ant_relay = 2;           % number of antenna of relay
num_phase_relay = 4;

M = 2;                     % Signal constellation size   
k = log2(M);                % No of bits per symbol

if M==2
   norfac =  1;
elseif M==4
   norfac =  1/sqrt(2);
elseif M==8
   norfac =  1/sqrt(6);
elseif M==18
    norfac =  1/sqrt(10);
elseif M==32
    norfac =  1/sqrt(26);
elseif M==64
    norfac =  1/sqrt(42);
elseif M==128
    norfac =  1/sqrt(106);
elseif M==256
    norfac =  1/sqrt(170);
end

%======= generate QAM constellation ======================
if M==2
   mapping = [-1 1]; 
elseif M==4
  mapping = [-1-j -1+j 1-j 1+j];
elseif M==8
  mapping = j*[-1-3*j -1-j -1+j -1+3*j 1-3*j 1-j 1+j 1+3*j];
elseif M==16
  mapping = [-3-3*j -3-j -3+j -3+3*j -1-3*j -1-j -1+j -1+3*j 1-3*j 1-j 1+j 1+3*j 3-3*j 3-j 3+j 3+3*j];
elseif M==32
  mapping = j*[-3-j*7 -3-j*5 -3-j*3 -3-j -3+j -3+3*j -3+5*j -3+7*j -1-j*7 -1-j*5 -1-j*3 -1-j -1+j -1+3*j -1+5*j -1+7*j 1-j*7 1-j*5 1-j*3 1-j 1+j 1+3*j 1+5*j 1+7*j 3-j*7 3-j*5 3-j*3 3-j 3+j 3+3*j 3+5*j 3+7*j];   
elseif M==64
  mapping = [-7-7*j -7-5*j -7-3*j -7-j -7+j -7+3*j -7+5*j -7+7*j -5-j*7 -5-j*5 -5-j*3 -5-j -5+j -5+3*j -5+5*j -5+7*j -3-j*7 -3-j*5 -3-j*3 -3-j -3+j -3+3*j -3+5*j -3+7*j -1-j*7 -1-j*5 -1-j*3 -1-j -1+j -1+3*j -1+5*j -1+7*j 1-j*7 1-j*5 1-j*3 1-j 1+j 1+3*j 1+5*j 1+7*j 3-j*7 3-j*5 3-j*3 3-j 3+j 3+3*j 3+5*j 3+7*j 5-j*7 5-j*5 5-j*3 5-j 5+j 5+3*j 5+5*j 5+7*j 7-j*7 7-j*5 7-j*3 7-j 7+j 7+3*j 7+5*j 7+7*j];
elseif M==128
  mapping = j*[-7-15*j -7-13*j -7-11*j -7-9*j -7-7*j -7-5*j -7-3*j -7-j -7+j -7+3*j -7+5*j -7+7*j -7+9*j -7+11*j -7+13*j -7+15*j -5-15*j -5-13*j -5-11*j -5-9*j -5-7*j -5-5*j -5-3*j -5-j -5+j -5+3*j -5+5*j -5+7*j -5+9*j -5+11*j -5+13*j -5+15*j -3-15*j -3-13*j -3-11*j -3-9*j -3-7*j -3-5*j -3-3*j -3-j -3+j -3+3*j -3+5*j -3+7*j -3+9*j -3+11*j -3+13*j -3+15*j -1-15*j -1-13*j -1-11*j -1-9*j -1-7*j -1-5*j -1-3*j -1-j -1+j -1+3*j -1+5*j -1+7*j -1+9*j -1+11*j -1+13*j -1+15*j 7-15*j 7-13*j 7-11*j 7-9*j 7-7*j 7-5*j 7-3*j 7-j 7+j 7+3*j 7+5*j 7+7*j 7+9*j 7+11*j 7+13*j 7+15*j 5-15*j 5-13*j 5-11*j 5-9*j 5-7*j 5-5*j 5-3*j 5-j 5+j 5+3*j 5+5*j 5+7*j 5+9*j 5+11*j 5+13*j 5+15*j 3-15*j 3-13*j 3-11*j 3-9*j 3-7*j 3-5*j 3-3*j 3-j 3+j 3+3*j 3+5*j 3+7*j 3+9*j 3+11*j 3+13*j 3+15*j 1-15*j 1-13*j 1-11*j 1-9*j 1-7*j 1-5*j 1-3*j 1-j 1+j 1+3*j 1+5*j 1+7*j 1+9*j 1+11*j 1+13*j 1+15*j];   
elseif M==256
  mapping = [-15-15*j -15-13*j -15-11*j -15-9*j -15-7*j -15-5*j -15-3*j -15-j -15+j -15+3*j -15+5*j -15+7*j -15+9*j -15+11*j -15+13*j -15+15*j -13-15*j -13-13*j -13-11*j -13-9*j -13-7*j -13-5*j -13-3*j -13-j -13+j -13+3*j -13+5*j -13+7*j -13+9*j -13+11*j -13+13*j -13+15*j -11-15*j -11-13*j -11-11*j -11-9*j -11-7*j -11-5*j -11-3*j -11-j -11+j -11+3*j -11+5*j -11+7*j -11+9*j -11+11*j -11+13*j -11+15*j -9-15*j -9-13*j -9-11*j -9-9*j -9-7*j -9-5*j -9-3*j -9-j -9+j -9+3*j -9+5*j -9+7*j -9+9*j -9+11*j -9+13*j -9+15*j -7-15*j -7-13*j -7-11*j -7-9*j -7-7*j -7-5*j -7-3*j -7-j -7+j -7+3*j -7+5*j -7+7*j -7+9*j -7+11*j -7+13*j -7+15*j -5-15*j -5-13*j -5-11*j -5-9*j -5-7*j -5-5*j -5-3*j -5-j -5+j -5+3*j -5+5*j -5+7*j -5+9*j -5+11*j -5+13*j -5+15*j -3-15*j -3-13*j -3-11*j -3-9*j -3-7*j -3-5*j -3-3*j -3-j -3+j -3+3*j -3+5*j -3+7*j -3+9*j -3+11*j -3+13*j -3+15*j -1-15*j -1-13*j -1-11*j -1-9*j -1-7*j -1-5*j -1-3*j -1-j -1+j -1+3*j -1+5*j -1+7*j -1+9*j -1+11*j -1+13*j -1+15*j 15-15*j 15-13*j 15-11*j 15-9*j 15-7*j 15-5*j 15-3*j 15-j 15+j 15+3*j 15+5*j 15+7*j 15+9*j 15+11*j 15+13*j 15+15*j 13-15*j 13-13*j 13-11*j 13-9*j 13-7*j 13-5*j 13-3*j 13-j 13+j 13+3*j 13+5*j 13+7*j 13+9*j 13+11*j 13+13*j 13+15*j 11-15*j 11-13*j 11-11*j 11-9*j 11-7*j 11-5*j 11-3*j 11-j 11+j 11+3*j 11+5*j 11+7*j 11+9*j 11+11*j 11+13*j 11+15*j 9-15*j 9-13*j 9-11*j 9-9*j 9-7*j 9-5*j 9-3*j 9-j 9+j 9+3*j 9+5*j 9+7*j 9+9*j 9+11*j 9+13*j 9+15*j 7-15*j 7-13*j 7-11*j 7-9*j 7-7*j 7-5*j 7-3*j 7-j 7+j 7+3*j 7+5*j 7+7*j 7+9*j 7+11*j 7+13*j 7+15*j 5-15*j 5-13*j 5-11*j 5-9*j 5-7*j 5-5*j 5-3*j 5-j 5+j 5+3*j 5+5*j 5+7*j 5+9*j 5+11*j 5+13*j 5+15*j 3-15*j 3-13*j 3-11*j 3-9*j 3-7*j 3-5*j 3-3*j 3-j 3+j 3+3*j 3+5*j 3+7*j 3+9*j 3+11*j 3+13*j 3+15*j 1-15*j 1-13*j 1-11*j 1-9*j 1-7*j 1-5*j 1-3*j 1-j 1+j 1+3*j 1+5*j 1+7*j 1+9*j 1+11*j 1+13*j 1+15*j];
end
%=========================================================
for i=1:2^k
    b=de2bi(i-1,k);     % b chua 2^k tu ma cua toan khong gian ma, moi tu ma dai k bit
    A(i)= mapping(1,bi2de(b(1:1:k),'left-msb')+1); %A la khong gian tin hieu
    A(i) = A(i);
end
%the codebook is ready
%===========================================================
% Generating data
datalen = 4;                 % number of variable in 2 phase
numsim = 1e+04;             % number of simulation bits
loop = round(numsim/datalen);%number of loop to achieve necesserary number of simulation bits
Am_node = sqrt(1/ant_node)*norfac;  % nomarlized transmit power factor/symbol duration at nodes = 1
Am_relay = sqrt(1/ant_relay)*norfac;  % nomarlized transmit power factor/symbol duration at relay = 1
% Eb/No definition
EbNodB = 0:2:22;              % Eb/No in dB
SNRdB = EbNodB + 10*log10(k);   % Change to Es/No
SNR = 10.^(SNRdB/10);
sigma = 1./sqrt(2*SNR); % noise deviation (No=sigma^2) per one dimension

%>>>>> SIMULATING BER <<<<<
num_exact=zeros(loop,length(SNRdB));
BER=zeros(loop,length(SNRdB));
for n=1:loop
    b = randint(datalen,1);	    % Tx binary data
    s = qammod(b,M);              % Modulate using M-PSK.
    %----------------
    clc;
    fprintf('Wait until loop =');disp(loop);
    fprintf('Current_loop =');disp(n);
    	% Generating quasi-static i.i.d. Rayleigh fading channels
    h = (randn(2,4)+j*randn(2,4))/sqrt(2);  
        
	for nn = 1:length(SNRdB) % Iteration over Eb/No
        %--------------
       z_node=sigma(nn)*(randn(2,4)+j*randn(2,4));% AWGN at nodes
       z_relay=sigma(nn)*(randn(1,4)+j*randn(1,4));% AWGN at relay
       %in phase 1: from nodes -> relay 
       H1 = Am_node*[h(1,1) h(1,2) h(2,1) h(2,2);
                     h(1,3) h(1,4) h(2,3) h(2,4)];
       X = [s(1);s(2);s(3);s(4)];
       noise1 = [z_relay(1,1);z_relay(1,2)];
       R1 = H1*X+noise1;
       
       %in phase 2: from nodes -> relay 
       H2 = Am_node*[conj(h(1,2)) conj(-h(1,1)) conj(h(2,2)) conj(-h(2,1));
                     conj(h(1,4)) conj(-h(1,3)) conj(h(2,4)) conj(-h(2,3))];
       X = [s(1); s(2);s(3);s(4)];
       noise2 = [z_relay(1,3);z_relay(1,4)];
       R2 = H2*X+noise2;
       
      % PNC in relay used ZF
       DD = [1 0 1 0;
            0 1 0 1;
            1 0 -1 0;
            0 1 0 -1];
       H= [H1;H2];
       R= [R1;R2];
       HH=H*inv(DD);
       G = pinv(HH);
       GG=G*G';
       Y = G*R;
       
          % PNC in relay used MMSE
%        DD = [1 0 1 0;
%             0 1 0 1;
%             1 0 -1 0;
%             0 1 0 -1];
%        H= [H1;H2];
%        R= [R1;R2];
%        HH=H*inv(DD);
%        G = (Am_node^2*HH'*HH+sigma(nn)^2*eye(4))^(-1)*HH'*Am_node^2;
%        GG=G*G';
%        Y = G*R;
          % Log Likehood Ratio - LLR decision method
       sigma1=GG(1,1)*2*sigma(nn)^2;
       sigma3=GG(3,3)*2*sigma(nn)^2;
       Y(1)=real(Y(1));Y(3)=real(Y(3));
       LLR1=exp(-Y(3)^2/(sigma3))*(exp(-(Y(1)-2)^2/(sigma1))+exp(-(Y(1)+2)^2/(sigma1)));
       LLR3=exp(-Y(1)^2/(sigma1))*(exp(-(Y(3)-2)^2/(sigma3))+exp(-(Y(3)+2)^2/(sigma3)));
       if LLR1>=LLR3
           y = -1;
       else y = 1;
       end
       %-----------
       sigma2=GG(2,2)*2*sigma(nn)^2;
       sigma4=GG(4,4)*2*sigma(nn)^2;
       Y(2)=real(Y(2));Y(4)=real(Y(4));
       LLR2=exp(-Y(4)^2/(sigma4))*(exp(-(Y(2)-2)^2/(sigma2))+exp(-(Y(2)+2)^2/(sigma2)));
       LLR4=exp(-Y(2)^2/(sigma2))*(exp(-(Y(4)-2)^2/(sigma4))+exp(-(Y(4)+2)^2/(sigma4)));
       if LLR2>=LLR4
           yy = -1;
       else yy = 1;
       end
       %-------------------% threshold decision
%        thr=1.0;
%        GG=G*G';
%        if GG(1,1)<GG(3,3)
%           y=-sign(abs(Y(1))-thr);
%        else y=-sign(thr-abs(Y(3)));
%        end
%        if GG(2,2)<GG(4,4)
%           yy=-sign(abs(Y(2))-thr);
%        else yy=-sign(thr-abs(Y(4)));
%        end
       %-------------------------

       %in phase 3,4: from relay -> nodes used Alamouti codes
       % for node1
       Hd = Am_relay*[h(1,1) h(1,3);
                      -conj(h(1,3)) conj(h(1,1));
                      h(1,2) h(1,4);
                      -conj(h(1,4)) conj(h(1,2))];
       Xd = [y; yy];
       noise_d = [z_node(1,1);z_node(1,2);z_node(1,3);z_node(1,4)];
       Rd = Hd*Xd+noise_d;
       D=Hd'*Rd;
       if real(D(1))>0
           x1=1;
       else x1=0;
       end
       if real(D(2))>0
           x2=1;
       else x2=0;
       end
       % De-XOR
       s3_hat = xor(b(1,1),x1);
       s4_hat = xor(b(2,1),x2);
       
       %-----------------
       % for node2
       Hd = Am_relay*[h(2,1) h(2,3);
                      -conj(h(2,3)) conj(h(2,1));
                      h(2,2) h(2,4);
                      -conj(h(2,4)) conj(h(2,2))];
       Xd = [y; yy];
       noise_d = [z_node(2,1);z_node(2,2);z_node(2,3);z_node(2,4)];
       Rd = Hd*Xd+noise_d;
       D=Hd'*Rd;
       if real(D(1))>0
           x1=1;
       else x1=0;
       end
       if real(D(2))>0
           x2=1;
       else x2=0;
       end
       % De-XOR
       s1_hat = xor(b(3,1),x1);
       s2_hat = xor(b(4,1),x2);
       %------------------------------
       % calculating error bit number
        if s1_hat==b(1,1);
            num_exact(n,nn)=num_exact(n,nn)+1;
        end
        if s2_hat==b(2,1);
            num_exact(n,nn)=num_exact(n,nn)+1;
        end
        if s3_hat==b(3,1);
            num_exact(n,nn)=num_exact(n,nn)+1;
        end
        if s4_hat==b(4,1);
            num_exact(n,nn)=num_exact(n,nn)+1;
        end
        num_error(n,nn)=datalen-num_exact(n,nn);
     %-----------
     BER(n,nn)=num_error(n,nn)/datalen;
    end
  %-----------     
  end
BER = sum(BER)/loop;
%save mat_PNC_MIMO_ShengliZhang_Alamouti_modif.mat BER;
semilogy(SNRdB,BER,'-rd');
legend ('ZF based MIMO-STBC-PNC (LLR)')
%legend ('ZF based MIMO-STBC-PNC (Selective)')
%legend ('MMSE based MIMO-STBC-PNC (LLR)')
%legend ('MMSE based MIMO-STBC-PNC (Selective)')
xlabel('Es / No [dB]');
ylabel('Bit Error Rate'); grid on;





