%===============================================================

clear all

% Write data to file
fileID1 = fopen('output_file_mod_zf.txt','r'); %ket qua tren modelsim
fileID2 = fopen('output_file_zf.txt','r'); %ketqua tren matlab
fileID3 = fopen('input_b1_b2.txt','r'); %gia tri b1b2

fileID4 = fopen('output_file_mod_mmse.txt','r'); %ket qua tren modelsim
fileID5 = fopen('output_file_mmse.txt','r'); %ketqua tren matlab

loop = 10000;
cnt_mod_zf = 0;
cnt_mat_zf = 0;

cnt_mod_mmse = 0;
cnt_mat_mmse = 0;

cnt = 0;
line = 0;
% Eb/No definition
datalen = 2; % Do dai tin hieu
k = 1;
EbNodB =0:2:26;                % Eb/No in dB
SNRdB = EbNodB + 10*log10(k);   % Change to Es/No
SNR = 10.^(SNRdB/10);
sigma = 1./sqrt(2*SNR);         % noise deviation (No=sigma^2) per one dimension

num_exact_mat_zf=zeros(loop,length(SNRdB)); %matlab chart
BER_mat_zf=zeros(loop,length(SNRdB));

num_exact_mat_mmse=zeros(loop,length(SNRdB)); %matlab chart
BER_mat_mmse=zeros(loop,length(SNRdB));
%----------------------------------------
num_exact_mod_zf=zeros(loop,length(SNRdB)); %modelsim chart
BER_mod_zf=zeros(loop,length(SNRdB));

num_exact_mod_mmse=zeros(loop,length(SNRdB)); %modelsim chart
BER_mod_mmse=zeros(loop,length(SNRdB));

for n=1:loop   

    %----------------
    %clc;
    %fprintf('Wait until loop =');disp(loop);
    %fprintf('Current_loop =');disp(n);
   	for nn = 1:length(SNRdB) % Iteration over Eb/No
        %--------------
        mod_zf = fscanf(fileID1,'%d',[1 2]);
        mat_zf = fscanf(fileID2,'%d',[1 2]);
        b1b2 = fscanf(fileID3,'%d',[1 2]);
        
        mod_mmse = fscanf(fileID4,'%d',[1 2]);
        mat_mmse = fscanf(fileID5,'%d',[1 2]);
%         % Moddelsim
%          fprintf('s1_mod = ');disp(mod(1, 1)); 
%          fprintf('s2_mod = ');disp(mod(1, 2));
%          fprintf('\n');
%         % Matlab
%         fprintf('s1_mat = ');disp(mat(1, 1)); 
%         fprintf('s2_mat = ');disp(mat(1, 2));
%         fprintf('\n');
    
%         % calculating error bit number
%         if (mat(1,1)>b1b2(1,1) || mat(1,1)<b1b2(1,1))
%              cnt_mat = cnt_mat + 1 ;
%         end
%     
%         if (mat(1,2)>b1b2(1,2) || mat(1,2)<b1b2(1,2))
%             cnt_mat = cnt_mat + 1 ;
%         end
%         %---------------------------------------------------
%         if (mod(1,1)>b1b2(1,1) || mod(1,1)<b1b2(1,1))
%              cnt_mod = cnt_mod + 1 ;
%         end
%     
%         if (mod(1,2)>b1b2(1,2) || mod(1,2)<b1b2(1,2))
%             cnt_mod = cnt_mod + 1 ;
%         end
%         %----------------------------------------------------
%         line = line + 1;
%         if (mat(1,1)>mod(1,1) || mat(1,1)<mod(1,1))
%              cnt = cnt + 1 ;
% %              fprintf('line= :');disp(line);
% %              fprintf('mat(1)= ');disp(mat(1,1));
% %              fprintf('mod(1)= \n');disp(mod(1,1));
%         end
%     
%         if (mat(1,2)>mod(1,2) || mat(1,2)<mod(1,2))
%             cnt = cnt + 1 ;
% %             fprintf('line= :');disp(line);
% %             fprintf('mat(2)= ');disp(mat(1,2));
% %             fprintf('mod(2)= \n');disp(mod(1,2));
%         end
%     if (mod(1,1)==mat(1,1));
%         cnt = cnt + 1 ;
%     end
        %Draw modelsim chart
         if mod_zf(1,1)==b1b2(1,1)
             num_exact_mod_zf(n,nn)=num_exact_mod_zf(n,nn)+1;
         end
         if mod_zf(1,2)==b1b2(1,2)
             num_exact_mod_zf(n,nn)=num_exact_mod_zf(n,nn)+1;
         end
        num_error_mod_zf=datalen-num_exact_mod_zf(n,nn);
     %-----------
     BER_mod_zf(n,nn)=num_error_mod_zf/datalen; 

     %-----------
     
         %Draw matlab chart
         if mat_zf(1,1)==b1b2(1,1)
             num_exact_mat_zf(n,nn)=num_exact_mat_zf(n,nn)+1;
         end
         if mat_zf(1,2)==b1b2(1,2)
             num_exact_mat_zf(n,nn)=num_exact_mat_zf(n,nn)+1;
         end
        num_error_mat_zf=datalen-num_exact_mat_zf(n,nn);
     %-----------
     BER_mat_zf(n,nn)=num_error_mat_zf/datalen; 
    
   
  %-----------
  
          %Draw modelsim chart
         if mod_mmse(1,1)==b1b2(1,1)
             num_exact_mod_mmse(n,nn)=num_exact_mod_mmse(n,nn)+1;
         end
         if mod_mmse(1,2)==b1b2(1,2)
             num_exact_mod_mmse(n,nn)=num_exact_mod_mmse(n,nn)+1;
         end
        num_error_mod_mmse=datalen-num_exact_mod_mmse(n,nn);
     %-----------
     BER_mod_mmse(n,nn)=num_error_mod_mmse/datalen; 

     %-----------
     
         %Draw matlab chart
         if mat_mmse(1,1)==b1b2(1,1)
             num_exact_mat_mmse(n,nn)=num_exact_mat_mmse(n,nn)+1;
         end
         if mat_mmse(1,2)==b1b2(1,2)
             num_exact_mat_mmse(n,nn)=num_exact_mat_mmse(n,nn)+1;
         end
        num_error_mat_mmse=datalen-num_exact_mat_mmse(n,nn);
     %-----------
     BER_mat_mmse(n,nn)=num_error_mat_mmse/datalen; 
    end
   
  %-----------     
end
  %Close file
  fclose(fileID1);
  fclose(fileID2);
  fclose(fileID3);
  fclose(fileID4);
  fclose(fileID5);

% fprintf('so bit loi matlab= \n');disp(cnt_mat);
% fprintf('so bit loi modelsim= \n');disp(cnt_mod); 
% fprintf('so bit khac nhau giua matlab va modelsim= \n');disp(cnt); 
BER_mod_zf = sum(BER_mod_zf)/loop;
BER_mat_zf = sum(BER_mat_zf)/loop;
BER_mod_mmse = sum(BER_mod_mmse)/loop;
BER_mat_mmse = sum(BER_mat_mmse)/loop;
%save mat_MIMO_SDM.mat BER;
semilogy(SNRdB,BER_mat_zf,'-rd', SNRdB,BER_mat_mmse,'-gs', SNRdB,BER_mod_zf,'-bd');% SNRdB,BER_mod_mmse,'-bs');
legend ('ZF based MIMO-SDM-PNC (LLR) Matlab', 'ZF new based MIMO-SDM-PNC (LLR) Modelsim', 'ZF based MIMO-SDM-PNC (LLR) Modelsim');% 'MMSE based MIMO-SDM-PNC (LLR) Modelsim');
xlabel('Es / No [dB]');
ylabel('Bit Error Rate'); grid on;




