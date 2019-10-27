function [Hhat_MMSE,C_MMSE,tau_p] = functionChannelEstimates_MMSE(H,R,...
    nbrOfSubcarriers,M,K,L,p)
%This function is used to generate the channel realizations and estimates 
%of these channels for all UEs in the entire network. The channels are 
%assumed to be correlated Rayleigh fading. The MMSE estimator is used.
%
%INPUT:
%
%H                 = M x nbrOfSubcarriers x K x L x L matrix with the
%                    channel realizations.
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices for all UEs in the network. R(:,:,k,j,l) is
%                    the correlation matrix for the channel between UE k 
%                    in cell j and the BS in cell l. This such matrix can
%                    either include the average channel gain or can be
%                    normalized arbitrarily.
%nbrOfSubcarriers  = Number of channel realizations (subcarriers)
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%
%Hhat_MMSE    = M x nbrOfSubcarriers x K x L x L matrix with the MMSE
%               channel estimates. The matrix Hhat_MMSE(:,n,k,j,l) is the
%               n:th channel estimate of the channel between UE k in cell j
%               and the BS in cell l.
%C_MMSE       = M x M x K x L x L matrix with estimation error correlation
%               matrices when using MMSE estimation. The matrix is
%               organized in the same way as R.
%tau_p        = Length of pilot sequences
%


%% Perform channel estimation


%Length of pilot sequences
tau_p = K;

%Store identity matrix of size M x M
eyeM = eye(M);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfSubcarriers,K,L) + 1i*randn(M,nbrOfSubcarriers,K,L));



%Prepare for MMSE estimation

%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M,nbrOfSubcarriers,K,L,L);

%Prepare to store estimation error correlation matrices
C_MMSE = zeros(M,M,K,L,L);


%% Go through all cells
for j = 1:L
    
    %Compute processed pilot signal for all UEs that use these pilots, according to (3.5)
    yp = sqrt(p)*tau_p*sum(H(:,:,:,:,j),4) + sqrt(tau_p)*Np(:,:,:,j);
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        
        %Go through all cells
        for l = 1:L
            
            %Check if the UE is active (inactive UEs have zero matrices)
            if trace(R(:,:,k,l,j))>0
                
                %Compute MMSE estimate of channel between BS l and UE k in
                %cell j using (3.9) in Theorem 3.1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = sqrt(p)*RPsi*yp(:,:,k);
                
                %Compute corresponding estimation error correlation matrix, using (3.11)
                C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*RPsi*R(:,:,k,l,j);
                
            end
            
        end
        
    end
    
end

