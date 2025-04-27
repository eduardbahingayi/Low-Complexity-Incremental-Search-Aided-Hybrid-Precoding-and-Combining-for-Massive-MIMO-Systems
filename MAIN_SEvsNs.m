clc;
clear all; %#ok<CLALL>
nChannel = 5e3; 
Nt = 64;   % TxAntennas
Nr = 64;    % RxAntennas
Ncl = 10;  % Number of clusters
Nray = 1;  % Number of rays
NRFmax = Ncl*Nray;
Ns = 4;     % Number of Data streams
NStreams = [Ns:NRFmax]; %Number of RF Chains
% Rate allocation vetors 
C1 = zeros(1,length(NStreams));
C2 = zeros(1,length(NStreams));
C3 = zeros(1,length(NStreams));
C4 = zeros(1,length(NStreams));
C5 = zeros(1,length(NStreams));
for SFind=1:length(NStreams)
    NRFr  = NStreams(SFind);% Number of RF Chains 
    NRFt = NRFr;
    SNR_dB = 0;  SNR = 10.^(SNR_dB/10)/Ns;
    %Output the progress
    %Output the progress
    disp(['Progress: snr = ' num2str(SFind) ' realizations.'])
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
        %% %%%%%%%%%%%%%%% proposed method  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [FRF]     = FastAntSelTX(ABS,SNR,H,NRFt);
        [WRF]     = FastAntSelRX(AMS,SNR,H*FRF,NRFt,NRFr);
        [FBB,WBB] = Baseband(WRF,H,FRF,Ns);                 
        F = FRF* FBB;
        W = WRF* WBB;
        temp5 = temp5 + log2(det(eye(Ns) + SNR * pinv(W) * H * (F * F') * H' * W));
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end        
    C1(SFind) = real(temp1/nChannel);
    C2(SFind) = real(temp2/nChannel);     
    C3(SFind) = real(temp3/nChannel);
    C4(SFind) = real(temp4/nChannel);     
    C5(SFind) = real(temp5/nChannel);  
end
%% figure;
hold on
plot( NStreams,C1,'r-p','Linewidth',1.5);
hold on
plot( NStreams,C2,'g->','Linewidth',1.5);
hold on
plot( NStreams,C3,'b--<','Linewidth',1.5);
hold on
plot( NStreams,C4,'Marker','d','Linewidth',1.5,'Color',[0 0.447058826684952 0.74117648601532]);
hold on
plot( NStreams,C5,'m-o','Linewidth',1.5);
hold off
legend({'Optimal unconstrained','PE-AltMin [15]','OMP [11]','CS-BHP [20]','Proposed algorithm'}...
    ,'Location','Northwest','FontSize',14,'FontWeight','normal','FontName' , 'Times New Roman');
xlabel('Number of RF chains (Mt=Mr)')
ylabel('Spectral efficiency (bps/Hz)')
grid on
box on
        