close all; 
clear all;
Nt=4;  % Number of TX
Nr=4;   %Number of Rx
Nchan=1e3;  %no of symbols
M = 4;
% d2b_bits = 2;
 
snrdb=0:25;%SNR in dB
snr = 10.^(0.1 .* snrdb);
%channel coeff
hsr= sqrt(1/2).*(((randn(Nr,Nt,Nchan))+1i.*(randn(Nr,Nt,Nchan))));
hrd= sqrt(1/2).*(((randn(Nr,Nt,Nchan))+1i.*(randn(Nr,Nt,Nchan))));
hsd= sqrt(1/2).*(((randn(Nr,Nt,Nchan))+1i.*(randn(Nr,Nt,Nchan))));
% random Tx antenna sequence %% [1, Nt] --> random integers btw 1->Nt 
Nb = log2(M*Nt);
x=randi([0,1],1,Nchan*Nb);  %% generating sequence for Antenna and QAM
x_1 = reshape(x,[], Nb);

antenna_bits = x_1(:,1:log2(Nt)); %% antenna id
q_bits = x_1(:,(log2(Nt)+1):Nb);  %% QAM data


antenna_dec = (bi2de(antenna_bits,'left-msb') +1); %% note: +1 is very important - antenna number start from 1 not 0

% QAM data generation
%b = randi([0 1],Nchan*log2(M),1);%Generate 0s and 1s
Qtx=qammod(q_bits',M,'bin','InputType','bit');                                                     %changed to transpose here



% AWGN
%n=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));
%n=zeros(Nr,Nchan);
for l=1:length(snrdb)
    l
%carrier signal
hxsr=zeros(Nr,Nchan);
hxrd=zeros(Nr,Nchan);
hxsd=zeros(Nr,Nchan);


% %%channel changes btw every transmission of bit
% for i=1:Nchan
%     
%     hxsr(:,i)=hsr(:,antenna_dec(i),i).*Qtx(i); %% : --> entire column, x(i) --> antenna , i --> symbol num
% end

%Transmission of signal with noise
%y=hx+n;

%Initializing loop for snr

    
n1=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));%getting noise inside loop
n2=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));
n3=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));
% n1=zeros(Nr,Nchan);
% n2=zeros(Nr,Nchan);
% n3=zeros(Nr,Nchan);
% Transmission with SNR data



%detection
%hx1=zeros(Nr,1,Nchan);



% generating QAM samples for comparision
QAM_init = 0:M-1;
Qsamp = qammod(QAM_init,M,'bin','InputType','integer');







%Source to Relay


for i=1:Nchan
    
    hxsr(:,i)=hsr(:,antenna_dec(i),i).*Qtx(i); %% : --> entire column, x(i) --> antenna , i --> symbol num
end

ysr=sqrt(snr(l)).*hxsr+n1;

QAM_data_sr = [];
ANT_data_sr = [];

for i = 1:Nchan  %% Iterating through all the symbols
   
    d = zeros(M,Nt);
    for j = 1:Nt    %% computing distance with all the possible data combo
        for k=1:M
            d(k,j) = norm(ysr(:,i) - sqrt(snr(l))*hsr(:,j,i)*Qsamp(k),'fro').^2; %% finding the distance of rx bit and actual bit from pilot sequence
        end
    end
    
    small_d= min(d(:)); %% estimating rx signal to data with min distance %% finds the index of min value
    [pos_qam, pos_ant] = find(small_d == d);
    temp = qamdemod(Qsamp(pos_qam),M,'bin','OutputType','bit');
    
    QAM_data_sr = [QAM_data_sr; temp'];
    ANT_data_sr = [ANT_data_sr pos_ant-1];
        
end



%relay to destination


Qtxrd=qammod(QAM_data_sr',M,'bin','InputType','bit');
%antenna_rd= (bi2de(antenna_bits,'left-msb') +1);




for i=1:Nchan
    
    hxrd(:,i)=hrd(:,ANT_data_sr(i)+1,i).*Qtxrd(i); %% : --> entire column, x(i) --> antenna , i --> symbol num
end

yrd=sqrt(snr(l)).*hxrd+n2;


QAM_data_rd = [];
ANT_data_rd = [];

for i = 1:Nchan  %% Iterating through all the symbols
   
    d = zeros(M,Nt);
    for j = 1:Nt    %% computing distance with all the possible data combo
        for k=1:M
            d(k,j) = norm(yrd(:,i) - sqrt(snr(l))*hrd(:,j,i)*Qsamp(k),'fro').^2; %% finding the distance of rx bit and actual bit from pilot sequence
        end
    end
    
    small_d= min(d(:)); %% estimating rx signal to data with min distance %% finds the index of min value
    [pos_qam, pos_ant] = find(small_d == d);
    temp = qamdemod(Qsamp(pos_qam),M,'bin','OutputType','bit');
    
    QAM_data_rd = [QAM_data_rd; temp'];
    ANT_data_rd = [ANT_data_rd pos_ant-1];
        
end


%Source to destination



QAM_data_sd = [];
ANT_data_sd = [];



%Qtxsd=qammod(QAM_data',M,'bin','InputType','bit');
%antenna_rd= (bi2de(antenna_bits,'left-msb') +1);




for i=1:Nchan
    
    hxsd(:,i)=hsd(:,antenna_dec(i),i).*Qtx(i); %% : --> entire column, x(i) --> antenna , i --> symbol num
end

ysd=sqrt(snr(l)).*hxsd+n3;



for i = 1:Nchan  %% Iterating through all the symbols
   
    d = zeros(M,Nt);
    for j = 1:Nt    %% computing distance with all the possible data combo
        for k=1:M
            d(k,j) = norm(ysd(:,i) - sqrt(snr(l))*hsd(:,j,i)*Qsamp(k),'fro').^2; %% finding the distance of rx bit and actual bit from pilot sequence
        end
    end
    
    small_d= min(d(:)); %% estimating rx signal to data with min distance %% finds the index of min value
    [pos_qam, pos_ant] = find(small_d == d);
    temp = qamdemod(Qsamp(pos_qam),M,'bin','OutputType','bit');
    
    QAM_data_sd = [QAM_data_sd; temp'];
    ANT_data_sd = [ANT_data_sd pos_ant-1];
        
end



%selection combining
for i = 1:Nchan
m=1;
        for j=1:Nt
            for a=1:Nt
                for q=1:M
                    for qcap=1:M
                if (j~=a)&&(q~=qcap)
                    gammasr(m)=norm(hsr(:,j,i).*Qsamp(q)-hsr(:,a,i).*Qsamp(qcap),'fro').^2;
                    m=m+1;
                end
                    end
                end
            end
        end
        gammasr1(i)=min(gammasr);

        
        
        m=1;
        for j=1:Nt
            for a=1:Nt
                for q=1:M
                    for qcap=1:M
                if (j~=a)&&(q~=qcap)
                    gammard(m)=norm(hrd(:,j,i).*Qsamp(q)-hrd(:,a,i).*Qsamp(qcap),'fro').^2;
                    m=m+1;
                end
                    end
                end
            end
        end
        gammard1(i)=min(gammard);
        
        
        
        
        m=1;
        for j=1:Nt
            for a=1:Nt
                for q=1:M
                    for qcap=1:M
                if (j~=a)&&(q~=qcap)
                    gammasd(m)=norm(hsd(:,j,i).*Qsamp(q)-hsd(:,a,i).*Qsamp(qcap),'fro').^2;
                    m=m+1;
                end
                    end
                end
            end
        end
        gammasd1(i)=min(gammasd);
        
        
        
        
        
end

for i = 1:Nchan
            if (gammasd1(i)>(min(gammasr1(i),gammard1(i))))
                ANT_data_final(i)=ANT_data_sd(i);
                QAM_data_final(i,:)=QAM_data_sd(i,:);
                %x_hatfinal(i)=x_hatsd(i);
            elseif(gammasd1(i)<(min(gammasr1(i),gammard1(i))))
                ANT_data_final(i)=ANT_data_rd(i);
                QAM_data_final(i,:)=QAM_data_rd(i,:);
                %x_hatfinal(i)=x_hatrd(i);
            end
end



ANT_bin_data_estimate_final = de2bi(ANT_data_final,'left-msb');
ANT_BER_final = size(find(ANT_bin_data_estimate_final - antenna_bits),1)/(log2(Nt*M)*Nchan);
  

QAM_BER_final = size(find(q_bits - QAM_data_final),1)/(log2(M*Nt)*Nchan);

BER_final(l) = log2(M)*ANT_BER_final + log2(Nt)*(QAM_BER_final);






ANT_bin_data_estimate_rd = de2bi(ANT_data_rd,'left-msb');
ANT_BER_rd = size(find(ANT_bin_data_estimate_rd - antenna_bits),1)/(log2(Nt*M)*Nchan);
  

QAM_BER_rd = size(find(q_bits - QAM_data_rd),1)/(log2(M*Nt)*Nchan);

BER_rd(l) = log2(M)*ANT_BER_rd + log2(Nt)*(QAM_BER_rd);






ANT_bin_data_estimate_sd = de2bi(ANT_data_sd,'left-msb');
ANT_BER_sd = size(find(ANT_bin_data_estimate_sd - antenna_bits),1)/(log2(Nt*M)*Nchan);
  

QAM_BER_sd = size(find(q_bits - QAM_data_sd),1)/(log2(M*Nt)*Nchan);

BER_sd(l) = log2(M)*ANT_BER_sd + log2(Nt)*(QAM_BER_sd);
end

semilogy(snrdb,BER_final);
hold on;
semilogy(snrdb,BER_rd);
semilogy(snrdb,BER_sd);
title('SM');
xlabel('Eb/N0, (dB)');
ylabel('BER');
legend('Using Cooperative Communication','Relay to Destination','Source to Destination');
%ylim([10^-5 0.2])
% grid on
% hold on;
% semilogy(snrdb,serrSD);
% semilogy(snrdb,serrRD);

