function [SE_DGOB] = functionComputeSE_DGOB(L,M,H,K,nbrOfSubcarriers,Pmax)
%This function is used to compute SE of FDD MIMO system using Digital Grid 
%Of Beams approach as described in [1]
%
%References:
%[1] Jose Flordelis et al. Massive MIMO Performance - TDD Versus FDD: 
%What Do Measurements Say? IEEE Transactions on Wireless Communications,
%Vol. 17, No. 4, April 2018
%
%INPUT:
%
%L = Number of BSs
%M = Number of BS array antennas
%H = M x nbrOfSubcarriers x K x L x L matrix with the exact channel 
%realizations
%K = number of BS UEs
%nbrOfSubcarriers = Number of channel realizations (subcarriers)
%Pmax = Maximum downlink transmit power per BS
%
%OUTPUT:
%
%SE_DGOB = K x L matrix of spectral efficiencies [b/s/Hz] calculated using
%Digital Grid Of Beams approach (DGOB)
%

%Number of BS array beams
Mb = M;

%Codebook matrix
C = dftmtx(Mb);

%Number of beams for UE to select and report
N = 30;



%Output simulation progress
disp('Iterating throught D-GOB algorithm');

%% Go through all channel realizations
nbrOfSubcarriers = 5; %FIXME to speed up debugging
%Prepare to store DGOB spectral efficiencies
SE_DGOB = zeros(K,L);
SEn = zeros(K,L,nbrOfSubcarriers); %per-subcarrier
for n = 1:nbrOfSubcarriers
    
    %Output simulation progress
    disp([num2str(n) ' channel realizations out of ' num2str(nbrOfSubcarriers)]);
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel realizations from all UEs of BS j to BS j
        Hnj = reshape(H(:,n,:,j,j),[M K]);
        Hhatnj = zeros(M,K);
        w = zeros(M,K);
        
        %Go through all UEs
        for k = 1:K
            
            g = C'*Hnj(:,k); %formula (8) from [1]
            
            [~, gs_index] = sort(g,'descend');
            Q = find(gs_index<=N); %select N strongest entries in g
            
            B = C(:,Q); %extracting relevant beams from C
            
            gu = B'*g; %formula (9) from [1] %#ok<NASGU>

            %Channel matrix quantization
            %Solve using CVX
            cvx_begin quiet
            
            variable v(M,1) complex;
            
            minimize norm(B'*v-gu) %formula (10) from [1]
            
            cvx_end
            
            %Use obtained v as Hhat (estimates) matrix
            Hhatnj(:,k) = v;
            
        end
        
        Z = pinv(Hhatnj'); %columns of the Moore-Penrose pseudoinverse
        
        %Prepare to store SINR and achievable sum-rate values
        SINR = zeros(1,K);
        Cdgob = zeros(1,K);
        
        for k = 1:K
            
            if norm(Z(:,k))>1e-10
                
                %ZF precoding
                w(:,k) = Z(:,k)/norm(Z(:,k));
                
                %Calculate interference term of (12) in [1]
                interf = 0;
                for i = 1:K
                    if i~=k
                        interf = interf + Hnj(:,i)'*w(:,i);
                    end
                end
                
                K_active = size(find(H(1,n,:,j,j)),1);
                
                SINR(k) = (Pmax/K_active*abs(Hnj(:,k)'*w(:,k)))^2 / ...
                    (1+Pmax/K_active*abs(interf)^2); %formula (12) in [1]
                
                Cdgob(k) = log2(1+SINR(k)); %formula (13) in [1]
                
                SEn(k,j,n) = Cdgob(k);
                
            end
            
        end
        
    end
end

%Go through all cells
for j = 1:L
    %Go through all UEs in cell j
    for k = 1:K
       SE_DGOB(k,j) = sum(SEn(k,j,:))/nbrOfSubcarriers; %average over 
       %subcarriers     
    end 
end


end





