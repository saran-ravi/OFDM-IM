close all; 
clear all;
Nt=4;
Nr=2;
Nchan=1e5;

hsr= sqrt(1/2).*(((randn(Nr,Nt,Nchan))+1i.*(randn(Nr,Nt,Nchan))));
hrd= sqrt(1/2).*(((randn(Nr,Nt,Nchan))+1i.*(randn(Nr,Nt,Nchan))));
hsd= sqrt(1/2).*(((randn(Nr,Nt,Nchan))+1i.*(randn(Nr,Nt,Nchan))));

x=randi([1,Nt],1,Nchan);
%n=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));

snrdb=0:20;
snr = 10.^(0.1 .* snrdb);
%count=0;
for k=1:length(snrdb)
    k
    %snr = 10.^(0.1 .* count)
     n1=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan)))); % noise wasmade common
     n2=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));
      n3=sqrt(1/2).*(((randn(Nr,Nchan))+1i.*(randn(Nr,Nchan))));
     %Source to Relay
    hxsr=zeros(Nr,Nchan);
    for i=1:Nchan
        hxsr(:,i)=hsr(:,x(i),i);
    end
    ysr=(sqrt(snr(k)).*hxsr)+n1;
    %ysr1=(sqrt(snr(k)).*hxsr)+n1;
    %detection
    
    x_hatsr = zeros(1,Nchan);
    for i = 1:Nchan
        dsr = zeros(1,Nt);
        for j = 1:Nt
            dsr(j) = norm(ysr(:,i) - sqrt(snr(k))*hsr(:,j,i),'fro').^2; %snr added to pilot data
        end
        [~,x_hatsr(1,i)]= min(dsr,[],2);
    end
    
    
    %Relay to destination
    
    
    hxrd=zeros(Nr,Nchan);
    for i=1:Nchan
        hxrd(:,i)=hrd(:,x_hatsr(i),i);
    end
    yrd=(sqrt(snr(k)).*hxrd)+n2;
    
    %detection
    
    x_hatrd = zeros(1,Nchan);
%     for i = 1:Nchan
%         drd = zeros(1,Nt);
%         for j = 1:Nt
%             drd(j) = norm(yrd(:,i) - sqrt(snr(k))*hrd(:,j,i),'fro').^2; %snr added to pilot data
%         end
%         [~,x_hatrd(1,i)]= min(drd,[],2);
%     end
    
    
    
    %source to destination
    
    
    
    hxsd=zeros(Nr,Nchan);
    for i=1:Nchan
        hxsd(:,i)=hsd(:,x(i),i);
    end
    ysd=(sqrt(snr(k)).*hxsd)+n3;
    
    %detection
    
    x_hatsd = zeros(1,Nchan);
    for i = 1:Nchan
        
        
        
        dsr = zeros(1,Nt);
        for j = 1:Nt
            dsr(j) = norm(ysr(:,i) - sqrt(snr(k))*hsr(:,j,i),'fro').^2; %snr added to pilot data
        end
        [dminsr(i),x_hatsr(1,i)]= min(dsr,[],2);
        
        
        
        
        drd = zeros(1,Nt);
        for j = 1:Nt
            drd(j) = norm(yrd(:,i) - sqrt(snr(k))*hrd(:,j,i),'fro').^2; %snr added to pilot data
        end
        [dminrd(i),x_hatrd(1,i)]= min(drd,[],2);
        
        
        
        
        dsd = zeros(1,Nt);
        for j = 1:Nt
            dsd(j) = norm(ysd(:,i) - sqrt(snr(k))*hsd(:,j,i),'fro').^2; %snr added to pilot data
        end
        [dminsd(i),x_hatsd(1,i)]= min(dsd,[],2);
        
        %gammasr=zeros(1,(Nt)^2-Nt);
        m=1;
        for j=1:Nt
            for a=1:Nt
                if j~=a
                    gammasr(m)=norm(hsr(:,j,i)-hsr(:,a,i),'fro').^2;
                    m=m+1;
                end
            end
        end
        gammasr1(i)=min(gammasr);
        
        
        
        %gammard=zeros(1,(Nt-1)^2);
        m=1;
        for j=1:Nt
            for a=1:Nt
                if j~=a
                    gammard(m)=norm(hrd(:,j,i)-hrd(:,a,i),'fro').^2;
                    m=m+1;
                end
            end
        end
        gammard1(i)=min(gammard);
        
        
        
        %gammasd=zeros(1,(Nt-1)^2);
        m=1;
        for j=1:Nt
            for a=1:Nt
                if j~=a
                    gammasd(m)=norm(hsd(:,j,i)-hsd(:,a,i),'fro').^2;
                    m=m+1;
                end
            end
        end
        gammasd1(i)=min(gammasd);
        
        
        
       
        
        
        %dcompsrd(i)=min(dminsr(i),dminrd(i));
%         if(dminsd(i)>(min(dminsr(i),dminrd(i))))
%             x_hatfinal(i)=x_hatsd(i);
%         elseif(dminsd(i)<(min(dminsr(i),dminrd(i))))
%             x_hatfinal(i)=x_hatrd(i);
%         end

            if (gammasd1(i)>(min(gammasr1(i),gammard1(i))))
                x_hatfinal(i)=x_hatsd(i);
            elseif(gammasd1(i)<(min(gammasr1(i),gammard1(i))))
                x_hatfinal(i)=x_hatrd(i);
            end
        
        
    end
    
%     for i=1:Nchan
%         x_hatcomp(i)=min(x_hatsr(i),x_hatrd(i));
%     end
%     
%     x_hatfinal = zeros(1,Nchan);
    
%     for i=1:Nchan
% %         x_hatcomp(i)=min(x_hatsr(i),x_hatrd(i));
%           if(x_hatcomp(i)>dsd(i))
%           x_hatfinal=
%     end
    
     serrRD(k)=size(find(x - x_hatrd),2)/Nchan;
     serrSD(k)=size(find(x - x_hatsd),2)/Nchan;
    serr(k)=size(find(x - x_hatfinal),2)/Nchan;
    %count=count+1;
end
%hold on;
semilogy(snrdb,serr);

title('SSK');
xlabel('Eb/N0, (dB)');
ylabel('Serr');
%ylim([10^-5 0.2])
grid on
hold on;
semilogy(snrdb,serrRD);
semilogy(snrdb,serrSD);
legend('Using Cooperative Communication','Relay to Destination','Source to Destination');
% title('SSK');
% xlabel('Eb/N0, (dB)');
% ylabel('Serr');
% grid on