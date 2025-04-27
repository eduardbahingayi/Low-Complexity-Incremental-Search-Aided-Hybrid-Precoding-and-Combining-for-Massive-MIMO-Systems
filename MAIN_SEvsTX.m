clc;
clear all; %#ok<CLALL>
nChannel = 1e2; 
NRFt = 4; 
NRFr = NRFt;
Ns = NRFr;   % Users = NRF Chains
Ncl = 5;   % number of rays
Nray = 1;   % number of rays
B = 8;
SNR_dB = 0;  SNR = 10.^(SNR_dB/10)/Ns;
TxRxAntenna = [32 64 128 256 384 512];
% Rate allocation vetors 
C1 = zeros(1,length(TxRxAntenna));
C2 = zeros(1,length(TxRxAntenna));
C3 = zeros(1,length(TxRxAntenna));
C4 = zeros(1,length(TxRxAntenna));
C5 = zeros(1,length(TxRxAntenna));
C6 = zeros(1,length(TxRxAntenna));
for Nind=1:length(TxRxAntenna)
    Nt = TxRxAntenna(Nind); % TxAntennas
    Nr = Nt;    % RxAntennas
    %Output the progress
    disp(['Progress: L = ' num2str(Nind) ' realizations.'])
    temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0; temp5 = 0; temp6 = 0;
    for m = 1:nChannel
        %Generate channel matrix for m:th realization
        %for model 1&2 to give the same results Nrays in model 1 
        %equals to Ncl in Model 2
        [H,Fopt,Wopt,ABS,AMS,D] = mmWaveChannelmodel(Nt, Nr, Ns, Ncl);
        %% %%%%%%%%%%%%%%% conventional method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [F,W] = SVD_precoding(H,Ns);
        temp1 = temp1 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W)); 
        %% %%%%%%%%%%%%%%% PE_AltMin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FRF,FBB] = PE_AltMin(Fopt,NRFt);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [WRF,WBB] = PE_AltMin(Wopt,NRFr);
        F = FRF* FBB;
        W = WRF* WBB;
        temp2 = temp2 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W));   
        %% %%%%%%%%%%%%%%% OMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FRF, FBB] = OMP(Fopt, ABS, NRFt);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        [WRF, WBB] = OMP(Wopt, AMS, NRFr);
        F = FRF* FBB;
        W = WRF* WBB;
        temp3 = temp3 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W));    
        %% %%%%%%%%%%%%%% CS-BHP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [AtDFT] = DFT_Codebook(Nt);% DFT dictionaries generator
        [FRF, FBB] = OBMP(Fopt, AtDFT, NRFt, Ns);
        F = FRF* FBB;
        WMMSE = (1/sqrt(SNR)*inv(Fopt'*(H'*H)*Fopt + 1/SNR*eye(Ns))*Fopt'*H')'; % rx, mmse
        W = WMMSE;
        temp4 = temp4 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W));    
        %% %%%%%%%%%%%%%%% Hybridly conneted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FRF, FBB ] = HybridlyConnectedTX(H,Nt,Ns);
        FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
        F = FRF* FBB;
        [WRF, WBB] = HybridlyConnectedRX(H,Nr,Ns);
        W = WRF* WBB;
        temp5 = temp5 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W)); 
        %% %%%%%%%%%%%%%%% proposed method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FRF]     = FastAntSelTX(ABS,SNR,H,NRFt);
        [WRF]     = FastAntSelRX(AMS,SNR,H*FRF,NRFt,NRFr);
        [FBB,WBB] = Baseband(WRF,H,FRF,Ns);                 
        F = FRF* FBB;
        W = WRF* WBB;
        temp6 = temp6 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W));
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    end
    C1(Nind) = real(temp1/nChannel);
    C2(Nind) = real(temp2/nChannel); 
    C3(Nind) = real(temp3/nChannel);
    C4(Nind) = real(temp4/nChannel);
    C5(Nind) = real(temp5/nChannel);
    C6(Nind) = real(temp6/nChannel);
 end
figure;
plot(TxRxAntenna, C1,'r-p','Linewidth',1.5);
hold on
plot(TxRxAntenna, C2,'g->','Linewidth',1.5);
hold on
plot(TxRxAntenna, C3,'b--<','Linewidth',1.5);
hold on
plot(TxRxAntenna, C4,'Marker','d','Linewidth',1.5,'Color',[0 0.447058826684952 0.74117648601532]);
hold on
plot(TxRxAntenna, C5,'Marker','s','Linewidth',1.5,'Color',[0.36 0.08 0.18]);
hold on
plot(TxRxAntenna, C6,'m-o','Linewidth',1.5);
hold off
legend({'Optimal unconstrained','PE-AltMin [15]','OMP [11]','CS-BHP [20]','SIC-HBF [18]','Proposed algorithm '}...
    ,'Location','Northwest','FontSize',14,'FontWeight','normal','FontName' , 'Times New Roman');
xlabel('Number of antennas (Nt = Nr) ')
ylabel('Spectral efficiency (bps/Hz)')
grid on
box on
        