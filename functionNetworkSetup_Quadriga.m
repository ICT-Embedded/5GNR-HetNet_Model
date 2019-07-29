function [H,Rest,activeUEs,pilotPattern,H_SC,Rest_SC,activeUEs_SC,...
    pilotPattern_SC,SCpositions,SCindex,Hbuilder_SC,SCindex_rnd,size_hex] = ...
    functionNetworkSetup_Quadriga(L,SCdrop,Kdrop,Kdrop_SC,B,B_SC,...
    noiseVariancedBm,noiseVariancedBm_SC,Kmax,Kmax_SC,f,M,M_SC,...
    polarizations,polarizations_SC,center_frequency,center_frequency_SC,...
    nbrOfSubs, nbrOfSubs_SC)
%This function generates the channel realizations between UEs at random
%locations and the BSs/SCs. BSs and SCs with cylindrical arrays and channel  
%are generated using QuaDRiGa from the Fraunhofer Heinrich Hertz Institute. 
%
%The 38.901 Urban Macrocell NLOS scenario is used for channel modeling of 
%macro cells.
%The 38.901 Urban Microcell LOS scenario is used for channel modeling of 
%small cells.
%
%QuaDRiGa needs to be installed separately (http://www.quadriga-channel-
%model.de) 
%and is delivered with a separate license. This function has been tested 
%using QuaDRiGa version 2.0.0-664.
%
%INPUT:
%L                = Number of BSs / macro cells
%SCdrop           = Number of SCs to be dropped within radius around BS
%Kdrop            = Number of UEs (including SCs) to be dropped within 
%                   radius around BS
%Kdrop_SC         = Number of UEs to be dropped within radius around SC
%B                = Bandwidth in Hz of macro layer
%B_SC             = Bandwidth in Hz of small cell layer
%noiseVariancedBm = Noise variance in dBm of macro layer
%noiseVariancedBm_SC = Noise variance in dBm of small cell layer
%Kmax             = Maximum number of UEs served by a BS
%Kmax_SC          = Maximum number of UEs served by a SC
%f                = Pilot reuse factor (macro cells), giving pilot 
%                   length Kmax*f
%M                = Number of BS antennas
%M_SC             = Number of SC antennas
%polarizations    = Select number of antenna polarizations (1 or 2) for 
%                   macro layer
%polarizations_SC = Select number of antenna polarizations (1 or 2) for 
%                   small cell layer
%center_frequency = Center frequency for macro layer
%center_frequency_SC = Center frequency for small cell layer
%nbrOfSubs = number of subcarriers, macro layer
%nbrOfSubs_SC = number of subcarriers, small cell layer
%
%OUTPUT:
%H         = M x nSubcarriers x K x L x L matrix with the channel 
%            realizations over nSubcarriers subcarriers at one time instance
%Rest      = M x M x K x L x L matrix with estimates of the spatial
%            correlation matrices for all UEs in the network.
%            Rest(:,:,k,j,l) is the correlation matrix for the channel
%            between UE k in cell j and the BS in cell l.
%activeUEs = Kmax x L with zeros and ones. activeUEs(k,l)==1 means that
%            pilot k is used by a UE in cell l
%pilotPattern = pilot pattern for BSs
%H_SC = M_SC x nSubcarriers_SC x K_SC x L_SC x L_SC matrix with the channel 
%            realizations over nSubcarriers subcarriers at one time
%            instance for small cell layer
%Rest_SC = M_SC x M_SC x K_SC x L_SC x L_SC matrix with estimates of the 
%            spatial correlation matrices for all UEs in small cell layer
%activeUEs_SC = Kmax_SC x L_SC with zeros and ones. activeUEs_SC(k,l)==1 
%            means that pilot k is used by a UE in small cell l
%pilotPattern_SC = pilot pattern for SCs
%SCpositions = Positions of small cells (macrocell UEs)
%SCindex = Indexes of UEs in macro cell to be 'backhauls' for small cells
%Hbuilder_SC = Kmax_SC x L_SC x L_SC matrix, with elements of custom 
%LinkGeometry class, which contains information about arrays, clusters and 
%angles/elevations of arrival/departure for the links
%SCindex_rnd = array of SC indexes corresponding to BS UEs which are SC 
%           backhauls 
%size_hex = lenght in [m] of macro hexagon's side
%

%% Create two new Quadriga layouts: one for macro layer (lay), one for 
%small cells (lay_SC)

%Output simulation progress
disp('Creating network layout with Quadriga');

%Set irrelevant parameters
s = qd_simulation_parameters;
s.sample_density = 1;
s.use_absolute_delays = 1;

%Set irrelevant parameters for SC
s_SC = qd_simulation_parameters;
s_SC.sample_density = 1;
s_SC.use_absolute_delays = 1;

%Set center frequency
%center_frequency = 4e9;
s.center_frequency = center_frequency;

%Set center frequency for SC
%center_frequency_SC = 28e9;
s_SC.center_frequency = center_frequency_SC;

%Number of subcarriers
nbrOfSubcarriers = nbrOfSubs;

%Number of subcarriers for SC
nbrOfSubcarriers_SC = nbrOfSubs_SC;

%Create new layout from general parameters
lay = qd_layout(s);

%Create new layout from general parameters for SC
lay_SC = qd_layout(s_SC);

%Generate BSs
lay.no_tx = L;

%Generate SCs
L_SC = SCdrop * L; % Total count of small cells in core macro cells
lay_SC.no_tx = L_SC;

%Set BS and SC height
BS_height = 25;
SC_height = 10;

%Set BS heights
lay.tx_position(3,:) = BS_height;

%Set SC heights
lay_SC.tx_position(3,:) = SC_height;


%% Deploy BSs

%Minimum distance between BSs and its UEs
minDistance = 35;

%Max and min distances between SC and its UE
maxDistance_SC = 50;
minDistance_SC = 10;

%Set the distance between BSs, in [m]
interSiteDistance = 200;

%Get the size of each cell hexagon (distance between the center and each 
%corner)
size_hex = interSiteDistance / sqrt(3);

%Horizontal and vertical spacing of hexagonal grid (between hexagon
%centers, i.e. base stations)
hs = sqrt(3) * size_hex;
vs = 3/2 * size_hex;

%Deploy BSs on hexagonal grid
% BSpositions = [...
%     3/2*hs + 1i*size_hex;... %1
%     5/2*hs + 1i*size_hex;... %2
%     7/2*hs + 1i*size_hex;... %3
%     hs + 1i*(size_hex+vs);... %4
%     2*hs + 1i*(size_hex+vs);... %5
%     3*hs + 1i*(size_hex+vs);... %6
%     4*hs + 1i*(size_hex+vs);... %7
%     1/2*hs  + 1i*(size_hex+2*vs) ;... %8
%     3/2*hs + 1i*(size_hex+2*vs);... %9
%     5/2*hs + 1i*(size_hex+2*vs);... %10
%     7/2*hs + 1i*(size_hex+2*vs);... %11
%     9/2*hs + 1i*(size_hex+2*vs);... %12
%     hs + 1i*(size_hex+3*vs);... %13
%     2*hs + 1i*(size_hex+3*vs);... %14
%     3*hs + 1i*(size_hex+3*vs);... %15
%     4*hs + 1i*(size_hex+3*vs);... %16
%     3/2*hs + 1i*(size_hex+4*vs);... %17
%     5/2*hs + 1i*(size_hex+4*vs);... %18
%     7/2*hs + 1i*(size_hex+4*vs)]; %19

%Small grid - for faster debugging
BSpositions = [...
    3/2*hs + 1i*size_hex;... %1
    5/2*hs + 1i*size_hex;... %2
    hs + 1i*(size_hex+vs);... %4
    2*hs + 1i*(size_hex+vs);... %5
    3*hs + 1i*(size_hex+vs);... %6
    3/2*hs + 1i*(size_hex+2*vs);... %9
    5/2*hs + 1i*(size_hex+2*vs)]; %10    

for j = 1:length(BSpositions)
    
    lay.tx_position(1:2,j) = [real(BSpositions(j)); imag(BSpositions(j))];
    
end


%% Create a circular antenna array for each BS
M_V = 5; %Number of vertical antennas
M_H = M/M_V; %Number of antennas on each horizontal circle

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

if polarizations == 1
    PolarizationIndicator = 1; %Single polarization (vertical) antennas
elseif polarizations == 2
    PolarizationIndicator = 3; %Dual +/-45deg polarized antennas
end

%Compute height of array
arrayHeight = (M_V-1)*antennaSpacing*3e8/center_frequency;

circumference = M_H*antennaSpacing*3e8/center_frequency;
radius = circumference/(2*pi);
delta_angle = 2*pi/M_H;

%Output simulation progress
disp('Creating macro cell arrays with Quadriga');

%Go through all BSs
for b = 1:L

    %Create rectangular array of size M_V x M_H
    lay.tx_array(b) = qd_arrayant('3gpp-3d', M_V, M_H, center_frequency,...
        PolarizationIndicator);

    %Place antennas on a circle and rotate radiation patters
    for i = 1:M_V
        for j = 1:M_H
            indices = (i-1)*M_H + j;
            angle = (j-1)*delta_angle;
            lay.tx_array(b).rotate_pattern(rad2deg(angle), 'z', indices, 0);
            lay.tx_array(b).element_position(1, indices) = radius*cos(angle);
            lay.tx_array(b).element_position(2, indices) = radius*sin(angle);
            lay.tx_array(b).element_position(3, indices) = ...
                (i-1)*antennaSpacing*3e8/center_frequency - arrayHeight/2;
        end
    end
end

%Plot BS cylindrical array
figure(41);
BS_txarray = phased.ConformalArray(...
    'ElementPosition',lay.tx_array(1).element_position,...
    'ElementNormal',zeros(2,M));
viewArray(BS_txarray)


%% Compute wrapping positions

%Compute seven (six + core) alternatives of the BS locations when using 
%wrap around for hexagonal grid
% wrapLocations = [...
%     0 + 1i*0 ...
%     7/2*hs - 1i*3*vs ...
%     4*hs + 1i*2*vs ...
%     1/2*hs + 1i*5*vs ...
%     -7/2*hs + 1i*3*vs ...
%     -4*hs - 1i*2*vs ...
%     -1/2*hs - 1i*5*vs];

%wrapping for small hex grid
wrapLocations = [...
    0 + 1i*0 ...
    2*hs - 1i*2*vs ...
    5/2*hs + 1i*1*vs ...
    1/2*hs + 1i*3*vs ...
    -2*hs + 1i*2*vs ...
    -5/2*hs - 1i*1*vs ...
    -1/2*hs - 1i*3*vs];


%Compute the exact dimension of the square/hexagon where the users are located
maxDistance = vs; %could be interSiteDistance

%Calculate BSpositionsWrapped
BSpositionsWrapped = zeros(L,length(wrapLocations));
for l = 1:L
    BSpositionsWrapped(l,:) = BSpositions(l) + wrapLocations;
end

%% Distribution UEs and SCs in macro cells

%Prepare to put out UEs in the cells
UEpositions = zeros(Kdrop,L);
UEpositionsWrapped = zeros(Kdrop,L,length(wrapLocations));
perBS = zeros(L,1);

%Select UEs which would be small cells
SCindex = zeros(SCdrop, L); % Indexes of UEs in macro cell to be 
%'backhauls' for small cells
SCpositions = zeros(SCdrop, L); % Positions of small cells (macrocell UEs)
SCpositionsWrapped = zeros(SCdrop,L,length(wrapLocations));

%Go through all the cells
for l = 1:L
    
    %Put out K UEs in the cell, uniformly at random. The procedure is
    %iterative since UEs that do not satisfy the minimum distance are
    %replaced with new UEs
    while perBS(l)<Kdrop
        
        %Put out users
        UEremaining = Kdrop-perBS(l);
        posX = rand(UEremaining,1)*maxDistance - maxDistance/2;
        posY = rand(UEremaining,1)*maxDistance - maxDistance/2;
        posXY = posX + 1i*posY;
        
        %Keep those that satisfy the minimum distance
        posXY = posXY(abs(posXY)>=minDistance);
        
        %Store new UEs
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + ...
            BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
        
    end
    
    min_dist_btw_SC = 50;
    dist_btw_SC = zeros(SCdrop-1,1);
    while min(dist_btw_SC) < min_dist_btw_SC        
        SCindex(:,l) = randperm(Kdrop,SCdrop);
        SCpositions(:,l) = UEpositions(SCindex(:,l),l);
        for i = 1:length(SCpositions(:,l))-1
            dist_btw_SC(i) = norm(SCpositions(i+1,l) - SCpositions(i,l));
        end
    end
    if SCdrop==1 %case of only one SC per macro and empty dist_btw_SC
        SCindex(:,l) = randperm(Kdrop,SCdrop);
        SCpositions(:,l) = UEpositions(SCindex(:,l),l);
    end
    
    %Create alternative UE positions using wrap around
    for k = 1:Kdrop
        
        UEpositionsWrapped(k,l,:) = UEpositions(k,l) + wrapLocations;
        
    end
    
    for sc = 1:SCdrop
        
        SCpositionsWrapped(sc,l,:) = SCpositions(sc,l) + wrapLocations;
        
    end
    
end

%generate tx_positions for Quadriga's small cell layer layout
for l = 1:L_SC

    l_sc = mod((l-1),SCdrop) + 1; %index of SC under BS (1..Kdrop_SC)
    l_bs = floor((l-1)/SCdrop) + 1; %index of BS (1..L)
    
    lay_SC.tx_position(1:2,l) = [real(SCpositions(l_sc,l_bs)); ...
        imag(SCpositions(l_sc,l_bs))];
    
end

%% Create a circular antenna array for each SC

M_V_SC = 8; %Number of vertical antennas of SC
M_H_SC = M_SC/M_V_SC; %Number of antennas on each horizontal circle of SC
P_V_SC = 1; %Number of vertical panels of SC
P_H_SC = 1; %Number of horizontal panels of SC

%Define the antenna and panel spacing (in number of wavelengths)
antennaSpacing_SC = 1/2; %Half wavelength distance
panelSpacing_SC = 2; %Double wavelength distance

if polarizations_SC == 1
    PolarizationIndicator_SC = 1; %Single polarization (vertical) antennas
elseif polarizations_SC == 2
    PolarizationIndicator_SC = 3; %Dual +/-45deg polarized antennas
end

%Compute height of array
arrayHeight_SC = (M_V_SC-1)*antennaSpacing_SC*3e8/center_frequency_SC;

circumference_SC = M_H_SC*antennaSpacing_SC*3e8/center_frequency_SC;
radius_SC = circumference_SC/(2*pi);
delta_angle_SC = 2*pi/M_H_SC;

%Output simulation progress
disp('Creating small cell arrays with Quadriga');

%Go through all SCs
for b = 1:L_SC %loop over SCs
    
    %Create rectangular array of size M_V x M_H
    lay_SC.tx_array(b) = qd_arrayant('3gpp-mmw', M_V_SC, M_H_SC, ...
        center_frequency_SC, PolarizationIndicator_SC, 0, ...
        antennaSpacing_SC, P_V_SC, P_H_SC, panelSpacing_SC, panelSpacing_SC);

    %Place antennas on a circle and rotate radiation patters
    for i = 1:M_V_SC
        for j = 1:M_H_SC
            indices_SC = (i-1)*M_H_SC + j;
            angle_SC = (j-1)*delta_angle_SC;
            lay_SC.tx_array(b).rotate_pattern(rad2deg(angle_SC), 'z', ...
                indices_SC, 0);
            lay_SC.tx_array(b).element_position(1, indices_SC) = ...
                radius_SC*cos(angle_SC);
            lay_SC.tx_array(b).element_position(2, indices_SC) = ...
                radius_SC*sin(angle_SC);
            lay_SC.tx_array(b).element_position(3, indices_SC) = ...
                (i-1)*antennaSpacing_SC*3e8/center_frequency_SC - ...
                arrayHeight_SC/2;
        end
    end
end

%Plot SC cylindrical array
figure(42);
SC_txarray = phased.ConformalArray(...
    'ElementPosition',lay_SC.tx_array(1).element_position,...
    'ElementNormal',zeros(2,M_SC));
viewArray(SC_txarray)


%% Distributing users in small cells
%Prepare to put out UEs in small cells
UEpositions_SC = zeros(Kdrop_SC,L_SC);
UEpositionsWrapped_SC = zeros(Kdrop_SC,L_SC,length(wrapLocations));
perSC = zeros(L_SC);

%Go through all small cells
%for sc = 1:SCdrop %loop over small cells of one macro cell
%for l = 1:L %loop over macro cells
for l = 1:L_SC % global index of SC (1..L_SC)
    l_sc = mod((l-1),SCdrop) + 1; % index of SC under BS (1..Kdrop_SC)
    l_bs = floor((l-1)/SCdrop) + 1; % index of BS (1..L)
    
    %Put out K_SC UEs in every small cell, uniformly at random. The procedure
    %is iterative since UEs that do not satisfy the minimum distance are
    %replaced with new UEs
    while perSC(l)<Kdrop_SC
        
        %Put out users
        UEremaining_SC = Kdrop_SC-perSC(l);
        posX_SC = rand(UEremaining_SC,1)*maxDistance_SC - maxDistance_SC/2;
        posY_SC = rand(UEremaining_SC,1)*maxDistance_SC - maxDistance_SC/2;
        posXY_SC = posX_SC + 1i*posY_SC;
        
        %Keep those that satisfy the minimum distance
        posXY_SC = posXY_SC(abs(posXY_SC)>=minDistance_SC);
        
        %Store new UEs
        UEpositions_SC(perSC(l)+1:perSC(l)+length(posXY_SC),l) = ...
            posXY_SC + SCpositions(l_sc,l_bs);
        perSC(l) = perSC(l)+length(posXY_SC);
        
    end
    
    %Create alternative UE positions using wrap around
    for k_sc = 1:Kdrop_SC
        
        UEpositionsWrapped_SC(k_sc,l,:) = UEpositions_SC(k_sc,l) + ...
            wrapLocations;
        
    end
end



%% Plot network layout (all cells + all users)

figure(31);

aH = axes;
scatter(aH,real(BSpositions), imag(BSpositions),15,'green');
hold on;

scatter(aH,real(BSpositionsWrapped(:)), imag(BSpositionsWrapped(:)));
aH.YDir = 'reverse';
%aH.XLim = [0 1800];
%aH.YLim = [0 1600];

voronoi(real([BSpositions' BSpositionsWrapped(:)']'), ...
    imag([BSpositions' BSpositionsWrapped(:)']'));
scatter(aH,real(UEpositions(:)), imag(UEpositions(:)),15,'green');
scatter(aH,real(UEpositionsWrapped(:)), imag(UEpositionsWrapped(:)));

scatter(aH,real(SCpositionsWrapped(:)), ...
    imag(SCpositionsWrapped(:)),15,'red');
scatter(aH,real(SCpositionsWrapped(:)), ...
    imag(SCpositionsWrapped(:)),800,'red');

scatter(aH,real(UEpositionsWrapped_SC(:)), ...
    imag(UEpositionsWrapped_SC(:)),15,'black');
hold off;


%% Configure macro UEs

Ktotal = Kdrop*L*length(wrapLocations); %Total number of macro UEs, 
%including SCs.

%Define macro UEs heights
UE_heights = 1.5*ones(Kdrop,L);
for i=1:L
    UE_heights(SCindex(:,i),i) = SC_height;
end
UE_heightsWrapped = repmat(UE_heights,[1 1 length(wrapLocations)]);

%Figure out indicies of Quadriga's macro layer layout receivers, which are 
%backhauls of small cells
K_backhaul = [];
for k = 1:SCdrop
    for l = 1:L
        for w = 1:length(wrapLocations)
            add = (Kdrop*(l-1)+SCindex(k,l))+Kdrop*L*(w-1);
            K_backhaul = [K_backhaul add];
        end
    end
end

%Generate macro UEs
lay.no_rx = Ktotal;

%Define macro UE antennas
lay.rx_array = qd_arrayant('omni');
%if we want to change backhaul antennas for SCs
% for k=1:Ktotal
%     if ~ismember(k,K_backhaul) %mobile subscribers
%         lay.rx_array(1,k) = qd_arrayant('omni');
%     else %fixed subscribers
%         lay.rx_array(1,k) = qd_arrayant('ula8');
%     end
% end


%% Simulate channels between macro BSs and their UEs (including SCs)

% Randomly distribute UEs
lay.rx_position =  [real(UEpositionsWrapped(:))'; ...
    imag(UEpositionsWrapped(:))'; UE_heightsWrapped(:)'];

%Define tracks for each UE, assuming a fixed UE location
for k=1:Ktotal
    %alternative way
    if ~ismember(k,K_backhaul) %mobile subscribers
        %lay.track(1,k) = qd_track('linear',50);
        lay.track(1,k).generate('linear',50);
        lay.track(1,k).name = ['Rx' num2str(k)];
        lay.track(1,k).scenario = '3GPP_38.901_UMa_NLOS';
        %lay.track(1,k).initial_position = lay.rx_position(:,k);
        %lay.track(1,k).positions = t.positions + t.initial_position;
        if mod(k,5)==0
            lay.track(1,k).set_speed(8.33); % 20% users 30km/h
        end
    else %fixed/'backhaul' subscribers
        %lay.track(1,k) = qd_track('linear',0,0);
        lay.track(1,k).generate('linear',0,0);
        lay.track(1,k).name = ['Rx' num2str(k)];
        lay.track(1,k).scenario = '3GPP_38.901_UMa_NLOS';
        %lay.track(1,k).initial_position = lay.rx_position(:,k);
        %lay.track(1,k).positions = t.positions + t.initial_position;
        lay.track(1,k).no_snapshots = 1;
    end
end

%Generate pilot patterns
if f == 1
    
    pilotPattern = ones(L,1);
    
elseif f == 2
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]); %Only works 
    %for 16 BSs
    
    %for L other than 16 do a simple alternation
    if L ~= 16
        pilotPattern = repmat([1; 2],floor(L/2),1);
        if mod(L,2)~=0
            pilotPattern = [pilotPattern; 1];
        end
    end
    
end


%Randomize pilot allocation in each cell
randOrder = zeros(Kmax*f,L);

for j = 1:L
    
    randOrder(1+(pilotPattern(j)-1)*Kmax:pilotPattern(j)*Kmax,j) = ...
        randperm(Kmax)+(pilotPattern(j)-1)*Kmax;
    
end


%Compute variance and standard deviation of the noise
noiseVar = 10^(noiseVariancedBm/10);
noiseStd = sqrt(noiseVar);

%Output simulation progress
disp('Generating channel realizations with Quadriga');

%Prepare to store channel realizations
H = zeros(M,nbrOfSubcarriers,Kmax,L,L);
Rest = zeros(M,M,Kmax,L,L);
perBS = zeros(L,1);
activeUEs = zeros(Kmax,L);

UEindex_rnd = zeros(Kdrop,L);
UEbestBS_rnd = zeros(Kdrop,L);

%% Go through all macro cells
for j = 1:L
    
    %Output simulation progress
    disp([num2str(j) ' cells generated out of ' num2str(L)]);
    
    %Go through all UEs of macro cell j
    for k = 1:Kdrop
        
        Huser = zeros(M,nbrOfSubcarriers,1,1,L);
        Ruser = zeros(M,M,1,1,L);
        
        %Extract the channels from UE k of macro cell j to *all* BSs
        for l = 1:L
            
            [~,minr] = min(abs(UEpositionsWrapped(k,j,:)-BSpositions(l)));
            
            userind = k+(j-1)*Kdrop+(minr-1)*Kdrop*L;
            
            [ h_channel, ~ ] = lay.get_channels_seg(l, userind);
            Hextract = h_channel.fr(B, nbrOfSubcarriers);  
            
            % indicies of Hextract are [Rx-Antenna , Tx-Antenna , 
            %Carrier-Index]
            Huser(:,:,1,1,l) = reshape(mean(Hextract(:,:,:,1),1),...
                [M nbrOfSubcarriers])/noiseStd;
            Ruser(:,:,1,1,l) = diag(mean(abs(Huser(:,:,1,1,l)).^2,2) / ...
                noiseVar);
            
        end
        
        %Determine which BS should serve the UE
        [~,bestBS] = max(mean(sum(abs(Huser(:,:,1,1,:)).^2,1),2));
        
        %Check if the selected BS has pilots available
        if perBS(bestBS)<Kmax
            
            %Add the UE to the cell of the selected BS
            perBS(bestBS) = perBS(bestBS) + 1;
            H(:,:,randOrder(perBS(bestBS)+ ...
                (pilotPattern(bestBS)-1)*Kmax,bestBS),bestBS,:) = Huser;
            Rest(:,:,randOrder(perBS(bestBS)+ ...
                (pilotPattern(bestBS)-1)*Kmax,bestBS),bestBS,:) = Ruser;
            activeUEs(randOrder(perBS(bestBS)+ ...
                (pilotPattern(bestBS)-1)*Kmax,bestBS),bestBS) = 1;
            
            UEindex_rnd(k,j) = randOrder(perBS(bestBS)+ ...
                (pilotPattern(bestBS)-1)*Kmax,bestBS);
            UEbestBS_rnd(k,j) = bestBS;
            
        end
        
    end
    
end

%Store indexes of BS UEs which are SC backhauls, and corresponding SC indexes
for l = 1:L
    for k = 1:SCdrop
        k_upd = UEindex_rnd(SCindex(k,l),l);
        l_upd = UEbestBS_rnd(SCindex(k,l),l);
        SCindex_rnd(k_upd,l_upd) = SCdrop*(l-1)+k;
    end
end


%% Configure SC UEs

Ktotal_SC = Kdrop_SC*SCdrop*L*length(wrapLocations); %Total number of SC UEs

%Define macro SC UEs heights
UE_heights_SC = 1.5*ones(Kdrop_SC,L_SC);
UE_heightsWrapped_SC = repmat(UE_heights_SC,[1 1 length(wrapLocations)]);

%Generate SC UEs
lay_SC.no_rx = Ktotal_SC;

%Define SC UE antennas
lay_SC.rx_array = qd_arrayant('omni');

%% Simulate channels between SCs and their UEs

% Randomly distribute UEs
lay_SC.rx_position =  [real(UEpositionsWrapped_SC(:))'; ...
    imag(UEpositionsWrapped_SC(:))'; UE_heightsWrapped_SC(:)'];

%Define tracks for each UE, assuming a fixed UE location
for k=1:Ktotal_SC
    lay_SC.track(k).generate('linear',0,0); %Define a linear track 
    %consisting of only one position
    %lay_SC.track(k).name = ['Rx' num2str(k)];
    lay_SC.track(k).scenario = '3GPP_38.901_UMi_LOS'; %Select the Urban 
    %Microcell LOS scenario
    lay_SC.track(k).no_snapshots = 1;

end
%Generate pilot patterns
%Do reuse 1 only since interference in mmW is expected to be negligible
pilotPattern_SC = ones(L_SC,1);

%Randomize pilot allocation in each cell
randOrder_SC = zeros(Kmax_SC,L_SC);

for j = 1:L_SC
    
    randOrder_SC(1+(pilotPattern_SC(j)-1)* ...
        Kmax_SC:pilotPattern_SC(j)*Kmax_SC,j) = ...
        randperm(Kmax_SC)+(pilotPattern_SC(j)-1)*Kmax_SC;
    
end

%Compute variance and standard deviation of the noise
%noiseVariancedBm_SC = -84;
noiseVar_SC = 10^(noiseVariancedBm_SC/10);
noiseStd_SC = sqrt(noiseVar_SC);


%Prepare to store channel realizations
H_SC = zeros(M_SC,nbrOfSubcarriers_SC,Kmax_SC,L_SC,L_SC);
%Hbuilder_SC = zeros(Kmax_SC,L_SC,L_SC);
Rest_SC = zeros(M_SC,M_SC,Kmax_SC,L_SC,L_SC);
perSC = zeros(L_SC,1);
activeUEs_SC = zeros(Kmax_SC,L_SC);


%% Go through all small cells
for j = 1:L_SC %global index of SC (1..L_SC)
    
    j_sc = mod((j-1),SCdrop) + 1; %index of SC under BS (1..SCdrop)
    j_bs = floor((j-1)/SCdrop) + 1; %index of BS (1..L)
    
    %Output simulation progress
    disp([num2str(j) ' small cells generated out of ' num2str(L_SC)]);
    
    %Go through all UEs of small cell j
    for k = 1:Kdrop_SC
        
        Huser_SC = zeros(M_SC,nbrOfSubcarriers_SC,1,1,L_SC);
        %Huser_builder_SC = zeros(L_SC,L_SC);
        Ruser_SC = zeros(M_SC,M_SC,1,1,L_SC);
        
        %Extract the channels from UE k of small cell j to *all* SCs
        for l = 1:L_SC %global index of SC (1..L_SC)
            
            l_sc = mod((l-1),SCdrop) + 1; %index of SC under BS (1..Kdrop_SC)
            l_bs = floor((l-1)/SCdrop) + 1; %index of BS (1..L)
            
            [~,minr_SC] = min(abs(UEpositionsWrapped_SC(k,j,:) - ...
                SCpositions(l_sc,l_bs)));
            
            userind_SC = k+(j-1)*Kdrop_SC+(minr_SC-1)*Kdrop_SC*L_SC;
            
            [ h_channel_SC, h_builder_SC ] = ...
                lay_SC.get_channels_seg(l, userind_SC);
            Hextract_SC = h_channel_SC.fr(B_SC, nbrOfSubcarriers_SC);
            
            % indicies of Hextract are [Rx-Antenna , Tx-Antenna , Carrier-Index]
            Huser_SC(:,:,1,1,l) = reshape(Hextract_SC, ...
                [M_SC nbrOfSubcarriers_SC])/noiseStd_SC;
            Huser_builder_SC(l) = LinkGeometry(...
                h_builder_SC(l).tx_array,...
                h_builder_SC(l).rx_array(userind_SC),...
                h_builder_SC(l).NumClusters,...
                h_builder_SC(l).AoD(userind_SC,:),...
                h_builder_SC(l).EoD(userind_SC,:),...
                h_builder_SC(l).AoA(userind_SC,:),...
                h_builder_SC(l).EoA(userind_SC,:));
            Ruser_SC(:,:,1,1,l) = ...
                diag(mean(abs(Huser_SC(:,:,1,1,l)).^2,2)/noiseVar_SC);
            
        end
        
        %Determine which SC should serve the UE
        [~,bestSC] = max(mean(sum(abs(Huser_SC(:,:,1,1,:)).^2,1),2));
        
        %Check if the selected BS has pilots available
        if perSC(bestSC)<Kmax_SC
            
            %Add the UE to the cell of the selected BS
            perSC(bestSC) = perSC(bestSC) + 1;
            H_SC(:,:,randOrder_SC(perSC(bestSC)+ ...
                (pilotPattern_SC(bestSC)-1)*Kmax_SC,bestSC),bestSC,:) = ...
                Huser_SC;
            Hbuilder_SC(randOrder_SC(perSC(bestSC)+ ...
                (pilotPattern_SC(bestSC)-1)*Kmax_SC,bestSC),bestSC,:) = ...
                Huser_builder_SC; %K x L_SC x L_SC
            Rest_SC(:,:,randOrder_SC(perSC(bestSC)+ ...
                (pilotPattern_SC(bestSC)-1)*Kmax_SC,bestSC),bestSC,:) = ...
                Ruser_SC;
            activeUEs_SC(randOrder_SC(perSC(bestSC)+ ...
                (pilotPattern_SC(bestSC)-1)*Kmax_SC,bestSC),bestSC) = 1;
            
        end
        
    end
    
end

%plot received power for macro cell tier (sub-6 GHz)
sample_dist = 5;
x_min = -500;
x_max = 1300;
y_min = -600;
y_max = 1200;
rx_height = 1.5;
tx_power = 30;
[map, x_coords, y_coords] = lay.power_map('3GPP_38.901_UMa_NLOS','quick',...
sample_dist,x_min,x_max,y_min,y_max,rx_height,tx_power);
sum_map = zeros(size(map{1},1));
for i=1:length(map)
    sum_map = sum_map + sum(map{i},4);
end
P_db = 10*log10(sum_map);

lay.visualize([],[],0);
hold on; imagesc(x_coords, y_coords, P_db); hold off;
axis([x_min,x_max,y_min,y_max]);
caxis(max(P_db(:)) + [-20 0] );
colmap = colormap;
colormap(colmap*0.5+0.5);
set(gca,'layer','top');
colorbar('south');
title('Received power [dBm]');


%plot received power for small cell tier (above-6 GHz)
sample_dist = 5;
x_min = -500;
x_max = 1300;
y_min = -600;
y_max = 1200;
rx_height = 1.5;
tx_power = 30;
[map, x_coords, y_coords] = lay_SC.power_map('3GPP_38.901_UMi_LOS','quick',...
sample_dist,x_min,x_max,y_min,y_max,rx_height,tx_power);
sum_map = zeros(size(map{1},1));
for i=1:length(map)
    sum_map = sum_map + sum(map{i},4);
end
P_db = 10*log10(sum_map);

lay_SC.visualize([],[],0);
hold on; imagesc(x_coords, y_coords, P_db); hold off;
axis([x_min,x_max,y_min,y_max]);
caxis(max(P_db(:)) + [-20 0] );
colmap = colormap;
colormap(colmap*0.5+0.5);
set(gca,'layer','top');
colorbar('south');
title('Received power [dBm]');
