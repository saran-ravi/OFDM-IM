clc
clear all
N = 4;
L=2;
Na = 1; 

N_symb = 1e4;

M = 4;
d = 0:M-1;
const = qammod(d,M,'UnitAveragePower', true);
nb_c = log2(M)*Na;

index_pattern = nchoosek([1:N],Na);
J = 2^floor(log2(nchoosek(N,Na)));
nb_si = log2(J);

index_pattern = index_pattern(1:J,:); % First J patterns are only taken since patterns are selected based on bits and J must be a power of two.

m = nb_si + nb_c;
Eb = (N + L)/m; 

snrdb=0:20;%SNR in dB
snr = 10.^(0.1 .* snrdb);

% Noise variance
N0T = Eb./snr; 
N0F = (Na/N)*N0T;



%DFT matrix calculation
W = fft(eye(N));


const_pattern = zeros(Na,M^Na);

for l = 1:Na
    m=1;
    j = 1;
    while m <=M^Na

        for n = 1:M^(l-1)
            const_pattern(l,m) =  const(j);
            m = m+1;
        end

        if mod(j+1,M) == 0
            j = M;
        else
            j = mod(j+1,M);
        end

    end
end


search_space = zeros(N,J*M^Na);

ind_pat_ite = 1;
i = 1;
while i <= (J)*M^Na



    for j = 1:M^Na
        search_space(index_pattern(ind_pat_ite,:),i) = const_pattern(:,j);
        i = i+1;
    end

    ind_pat_ite = ind_pat_ite + 1;

end
% search_space
% size(search_space)

t = nchoosek(J*M^Na,2);


bits = randi([0 1],nb_si+nb_c, N_symb);
d_bits = bi2de(bits.','left-msb').';


for k=1:length(snrdb)
    k
d_bits_cap_sd = zeros(1,N_symb);
d_bits_cap_sr = zeros(1,N_symb);
d_bits_cap_rd = zeros(1,N_symb);
    for i=1:N_symb
%-------------------Source to destination

        %% Tx
        
        xf_sd= search_space(:,d_bits(i)+1); % xf in paper

        s_tilde_sd = (1/sqrt(Na)).*(W')*xf_sd; %eqn.(6) in paper
        s_prime_sd = [s_tilde_sd(end - L + 2:end,:) ; s_tilde_sd]; % adding cyclic prefix
        
        %% Channel
        gsd = sqrt(1/(2*L)).*(randn(L,1) + 1i.*randn(L,1));


        htsd= [gsd; zeros(N-L,1)];
        

        hfsd = W*htsd; % eqn.(16) in paper
        
        %noise
        n1t=sqrt(N0T(k)/2).*(randn(N,1) + 1i.*randn(N,1));
        
        n1f = (sqrt(Na)/N).*W*n1t; % noise passed through (sqrt(Na)/N).*W (fft block) at rx
        
        yfsd = diag(xf_sd)*hfsd + n1f;
        
        
        
        %% Rx
        eucld_d_sd = zeros(1,J*M^Na);
        for l = 1:J*M^Na
            temp_sd = diag(search_space(:,l))*hfsd;
            eucld_d_sd(l) = norm(yfsd - temp_sd,'fro')^2;
        end

        [~,b_sd]=min(eucld_d_sd,[],2);
        d_bits_cap_sd(i) = b_sd-1;
        
        
        %----------------Source to relay
        

 %% Tx
        
        xf_sr= search_space(:,d_bits(i)+1); % xf in paper

        s_tilde_sr = (1/sqrt(Na)).*(W')*xf_sr; %eqn.(6) in paper
        s_prime_sr = [s_tilde_sr(end - L + 2:end,:) ; s_tilde_sr]; % adding cyclic prefix
        
        %% Channel
        gsr = sqrt(1/(2*L)).*(randn(L,1) + 1i.*randn(L,1));


        htsr= [gsr; zeros(N-L,1)];
        

        hfsr = W*htsr; % eqn.(16) in paper
        
        %noise
        n2t=sqrt(N0T(k)/2).*(randn(N,1) + 1i.*randn(N,1));
        
        n2f = (sqrt(Na)/N).*W*n2t; % noise passed through (sqrt(Na)/N).*W (fft block) at rx
        
        yfsr = diag(xf_sr)*hfsr + n2f;
        
        
        
        %% Rx
        eucld_d_sr = zeros(1,J*M^Na);
        for l = 1:J*M^Na
            temp_sr = diag(search_space(:,l))*hfsr;
            eucld_d_sr(l) = norm(yfsr - temp_sr,'fro')^2;
        end

        [~,b_sr]=min(eucld_d_sr,[],2);
        d_bits_cap_sr(i) = b_sr-1;
        
        
        
        
        %----------------------------Relay to destination
        
         %% Tx
        
        xf_rd= search_space(:,d_bits_cap_sr(i)+1); % xf in paper

        s_tilde_rd = (1/sqrt(Na)).*(W')*xf_rd; %eqn.(6) in paper
        s_prime_rd = [s_tilde_rd(end - L + 2:end,:) ; s_tilde_rd]; % adding cyclic prefix
        
        %% Channel
        grd = sqrt(1/(2*L)).*(randn(L,1) + 1i.*randn(L,1));


        htrd= [grd; zeros(N-L,1)];
        

        hfrd = W*htrd; % eqn.(16) in paper
        
        %noise
        n3t=sqrt(N0T(k)/2).*(randn(N,1) + 1i.*randn(N,1));
        
        n3f = (sqrt(Na)/N).*W*n3t; % noise passed through (sqrt(Na)/N).*W (fft block) at rx
        
        yfrd = diag(xf_rd)*hfrd + n3f;
        
        
        
        %% Rx
        eucld_d_rd = zeros(1,J*M^Na);
        for l = 1:J*M^Na
            temp_rd = diag(search_space(:,l))*hfrd;
            eucld_d_rd(l) = norm(yfrd - temp_rd,'fro')^2;
        end

        [~,b_rd]=min(eucld_d_rd,[],2);
        d_bits_cap_rd(i) = b_rd-1;

        
        
        d_sd = zeros(1,t);
    d_sr = zeros(1,t);
    d_rd = zeros(1,t);

    o=1;
    for m = 1:J*M^Na
        for j = m+1:J*M^Na

            X = diag(search_space(:,m));
            X_cap = diag(search_space(:,j));

            d_sd(o) = norm((X - X_cap)*hfsd,'fro')^2; % delta given in paper below eqn.(21)
            d_sr(o) = norm((X - X_cap)*hfsr,'fro')^2; % delta given in paper below eqn.(21)
            d_rd(o) = norm((X - X_cap)*hfrd,'fro')^2; % delta given in paper below eqn.(21)
            
            o = o+1;
            
        end
    end
    
    delta_min_sd(i) = min(d_sd,[],2);
    delta_min_sr(i) = min(d_sr,[],2);
    delta_min_rd(i) = min(d_rd,[],2);
    
    if delta_min_sd(i) > min(delta_min_sr(i),delta_min_rd(i))
        
        %choose sd link
        d_bits_cap_sc(i)=d_bits_cap_sd(i);
    else
        %choose rd link
        d_bits_cap_sc(i)=d_bits_cap_rd(i);
    end


        



    end
    [d_bits; d_bits_cap_sd];
    %-------------------Using Selection Combining
        bits_re=reshape(bits,1,(nb_si+nb_c)*N_symb);
    
    for i=1:N_symb
    bits_rec_sc(:,i)=de2bi(d_bits_cap_sc(i),nb_si+nb_c,'left-msb');
    end
    
    bits_rec_re_sc=reshape(bits_rec_sc,1,(nb_si+nb_c)*N_symb);
    BER_sc(k)=size(find(bits_re-bits_rec_re_sc),2)/((nb_si+nb_c)*N_symb);
    %--------------S to D

    
    for i=1:N_symb
    bits_rec_sd(:,i)=de2bi(d_bits_cap_sd(i),nb_si+nb_c,'left-msb');
    end
    
    bits_rec_re_sd=reshape(bits_rec_sd,1,(nb_si+nb_c)*N_symb);
    BER_sd(k)=size(find(bits_re-bits_rec_re_sd),2)/((nb_si+nb_c)*N_symb);
    
    %-----------------------R to D
    
    for i=1:N_symb
    bits_rec_rd(:,i)=de2bi(d_bits_cap_rd(i),nb_si+nb_c,'left-msb');
    end
    
    bits_rec_re_rd=reshape(bits_rec_rd,1,(nb_si+nb_c)*N_symb);
    BER_rd(k)=size(find(bits_re-bits_rec_re_rd),2)/((nb_si+nb_c)*N_symb);
    
end

semilogy(snrdb,BER_sc);
xlabel('Eb/N0, (dB)');
ylabel('BER');
ylim([1e-3 1]);
hold on;
semilogy(snrdb,BER_sd);
semilogy(snrdb,BER_rd);
legend('Using Selection Combining','Source to Destination','Relay to Destination');






