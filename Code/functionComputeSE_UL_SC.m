function [SE_UL_SC] = functionComputeSE_UL_SC(H,Hhat,nbrOfSubcarriers,...
    M,K,L,Hbuilder_SC,Pmax_SC,center_frequency_SC)
%
%This function is used to compute UL spectral efficiency for small cells
%
%INPUT:
%
%H                 = M x nbrOfSubcarriers x K x L x L matrix with the
%                    exact channel realizations
%Hhat              = M x nbrOfSubcarriers x K x L x L matrix with the MMSE
%                    channel estimates
%nbrOfSubcarriers  = Number of channel realizations (subcarriers)
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%Hbuilder_SC       = Kmax_SC x L_SC x L_SC matrix, with elements of custom 
%LinkGeometry class, which contains information about arrays, clusters and 
%angles/elevations of arrival/departure for the links
%Pmax_SC           = SC transmit power
%center_frequency_SC = center frequency of small cell layer
%
%OUTPUT:
%
%SE_UL_SC = K x L matrix with spectral efficiencies [b/s/Hz] calculated 
%using OMP approach
%


%Prepare to store spectral efficiencies
SE_UL_SC = zeros(K,L);

%% Go through all channel realizations
for n = 1:nbrOfSubcarriers
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);

        %Go through all UEs in cell j
        for k = 1:K
            
            if norm(V_MR(:,k))>0
  
                SE_UL_SC(k,j) = SE_UL_SC(k,j) + ...
                    functionComputeSE_OMP(Hbuilder_SC(k,j,j),...
                    H(:,n,k,j,j), Pmax_SC/2, center_frequency_SC);
                
            end
            
        end
        
    end
    
end



SE_UL_SC = SE_UL_SC / nbrOfSubcarriers; %average over subcarriers
