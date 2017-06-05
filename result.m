%===============================================================

clear all

% Write data to file
fileID1 = fopen('output_file_mod.txt','r'); %ket qua tren modelsim
fileID2 = fopen('output_file.txt','r'); %ketqua tren matlab
fileID3 = fopen('input_b1_b2.txt','r'); %gia tri b1b2

loop = 1000;
cnt_mod = 0;
cnt_mat = 0;
cnt = 0;
line = 0;
% Eb/No definition
datalen = 2; % Do dai tin hieu
k = 1;
EbNodB =0:2:26;                % Eb/No in dB
SNRdB = EbNodB + 10*log10(k);   % Change to Es/No
SNR = 10.^(SNRdB/10);
sigma = 1./sqrt(2*SNR);         % noise deviation (No=sigma^2) per one dimension

num_exact_mat=zeros(loop,length(SNRdB)); %matlab chart
BER_mat=zeros(loop,length(SNRdB));
%----------------------------------------
num_exact_mod=zeros(loop,length(SNRdB)); %modelsim chart
BER_mod=zeros(loop,length(SNRdB));

for n=1:loop   

    %----------------
    %clc;
    %fprintf('Wait until loop =');disp(loop);
    %fprintf('Current_loop =');disp(n);
   	for nn = 1:length(SNRdB) % Iteration over Eb/No
        %--------------
        mod = fscanf(fileID1,'%d',[1 2]);
        mat = fscanf(fileID2,'%d',[1 2]);
        b1b2 = fscanf(fileID3,'%d',[1 2]);
%         % Moddelsim
%          fprintf('s1_mod = ');disp(mod(1, 1)); 
%          fprintf('s2_mod = ');disp(mod(1, 2));
%          fprintf('\n');
%         % Matlab
%         fprintf('s1_mat = ');disp(mat(1, 1)); 
%         fprintf('s2_mat = ');disp(mat(1, 2));
%         fprintf('\n');
    
        % calculating error bit number
        if (mat(1,1)>b1b2(1,1) || mat(1,1)<b1b2(1,1))
             cnt_mat = cnt_mat + 1 ;
        end
    
        if (mat(1,2)>b1b2(1,2) || mat(1,2)<b1b2(1,2))
            cnt_mat = cnt_mat + 1 ;
        end
        %---------------------------------------------------
        if (mod(1,1)>b1b2(1,1) || mod(1,1)<b1b2(1,1))
             cnt_mod = cnt_mod + 1 ;
        end
    
        if (mod(1,2)>b1b2(1,2) || mod(1,2)<b1b2(1,2))
            cnt_mod = cnt_mod + 1 ;
        end
        %----------------------------------------------------
        line = line + 1;
        if (mat(1,1)>mod(1,1) || mat(1,1)<mod(1,1))% || mat(1,2)>mod(1,2) || mat(1,2)<mod(1,2))
             cnt = cnt + 1 ;
%              fprintf('line= :');disp(line);
%              fprintf('mat(1)= ');disp(mat(1,1));
%              fprintf('mod(1)= \n');disp(mod(1,1));
        end
    
        if (mat(1,2)>mod(1,2) || mat(1,2)<mod(1,2))
            cnt = cnt + 1 ;
%             fprintf('line= :');disp(line);
%             fprintf('mat(2)= ');disp(mat(1,2));
%             fprintf('mod(2)= \n');disp(mod(1,2));
        end
%     if (mod(1,1)==mat(1,1));
%         cnt = cnt + 1 ;
%     end
        %Draw modelsim chart
         if mod(1,1)==b1b2(1,1)
             num_exact_mod(n,nn)=num_exact_mod(n,nn)+1;
         end
         if mod(1,2)==b1b2(1,2)
             num_exact_mod(n,nn)=num_exact_mod(n,nn)+1;
         end
        num_error_mod=datalen-num_exact_mod(n,nn);
     %-----------
     BER_mod(n,nn)=num_error_mod/datalen; 

     %-----------
     
         %Draw matlab chart
         if mat(1,1)==b1b2(1,1)
             num_exact_mat(n,nn)=num_exact_mat(n,nn)+1;
         end
         if mat(1,2)==b1b2(1,2)
             num_exact_mat(n,nn)=num_exact_mat(n,nn)+1;
         end
        num_error_mat=datalen-num_exact_mat(n,nn);
     %-----------
     BER_mat(n,nn)=num_error_mat/datalen; 
    end
   
  %-----------     
end
  %Close file
  fclose(fileID1);
  fclose(fileID2);
  fclose(fileID3);

fprintf('so bit loi matlab= \n');disp(cnt_mat);
fprintf('so bit loi modelsim= \n');disp(cnt_mod); 
fprintf('so bit khac nhau giua matlab va modelsim= \n');disp(cnt); 
BER_mod = sum(BER_mod)/loop;
BER_mat = sum(BER_mat)/loop;
%save mat_MIMO_SDM.mat BER;
semilogy(SNRdB,BER_mat,'-rd', SNRdB,BER_mod,'-sb');
legend ('MMSE based MIMO-SDM-PNC (LLR) Matlab', 'MMSE based MIMO-SDM-PNC (LLR) Modelsim');
xlabel('Es / No [dB]');
ylabel('Bit Error Rate'); grid on;




