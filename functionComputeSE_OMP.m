function [SE_OMP] = functionComputeSE_OMP(Hbuilder_SC, Hh, ...
    Pmax_SC, center_frequency_SC)
%
%This function is used to compute spectral efficiency of TDD MIMO system 
%with hybrid or analog beamforming in place. It uses the approach proposed 
%in [1] which decouples the optimizations for the precoding and combining 
%weights. It first uses the orthogonal matching pursuit algorithm to derive 
%the precoding weights. Once the precoding weights are computed, 
%the result is then used to obtain the corresponding combining weights.
%
%References:
%[1] Oma El Ayach, et al. Spatially Sparse Precoding in Millimeter wave 
%MIMO Systems, IEEE Transactions on Wireless Communications, Vol. 13, 
%No. 3, March 2014.
%[2] https://www.mathworks.com/help/phased/examples/introduction-to-hybrid-
%beamforming.html
%
%INPUT:
%
%Hbuilder_SC = Kmax_SC x L_SC x L_SC matrix, with elements of custom 
%LinkGeometry class, which contains information about arrays, clusters and 
%angles/elevations of arrival/departure for the links
%Hh = channel realization for the particular channel between UE k in cell j 
%and the BS in cell j for subcarrier n (H(:,n,k,j,j))
%Pmax_SC = SC transmit power
%center_frequency_SC = center frequency of small cell layer
%
%OUTPUT:
%
%SE_OMP = spectral efficiency for particular link, calculated using 
%Orthogonal Matching Pursuit approach (OMP)
%

Nt = max(size(Hh));
NtRF = 1;

Nr = min(size(Hh));
NrRF = 1;

Ns = 1;

fc = center_frequency_SC;

At = steervec(Hbuilder_SC.tx_array.element_position,rad2deg([Hbuilder_SC.AoD;...
    Hbuilder_SC.EoD]));
Ar = steervec(Hbuilder_SC.rx_array.element_position,rad2deg([Hbuilder_SC.AoA;...
    Hbuilder_SC.EoA]));

[Fbb,Frf] = helperOMPHybridPrecodingWeights(Hh,NtRF,Ns,At);
[Fbb,Frf,Wbb,Wrf] = helperOMPHybridWeights(Hh,NtRF,NrRF,Ns,At,Ar,1/(db2pow(sqrt(Pmax_SC))/2));

SE_OMP = helperComputeSpectralEfficiency(Hh,Fbb*Frf,Wrf*Wbb,Ns,db2pow(sqrt(Pmax_SC))/2);

%%Plot array patterns

% txarray1 = phased.ConformalArray(...
%     'ElementPosition',Hbuilder_SC.tx_array.element_position,...
%     'ElementNormal',ones(2,Nt));
% 
% txarray2 = phased.PartitionedArray(...
%     'Array',txarray1,...
%     'SubarraySelection',ones(NtRF,Nt),...
%     'SubarraySteering','Custom');

%(Conformal array) The beam pattern of the hybrid weights is shown below
% figure(72);
% pattern(txarray1,fc,-180:180,-90:90,...
%     'Type','powerdb',...
%     'Weights',Frf'*Fbb');

%(PartitionedArray) The beam pattern of the hybrid weights is shown below
% figure(73);
% pattern(txarray2,fc,-180:180,-90:90,...
%     'Type','efield',...
%     'ElementWeights',Frf'*Fbb');