function [H,Fopt,Wopt,ABS,AMS,d]= mmWaveChannelmodel(Nt, Nr, Ns, Ncl)
fc = 28e9; % Frequencey 
lamada = 3e8/fc; % wavelegenth;
gamma = sqrt(Nr*Nt/Ncl);
alpha = (randn(1,Ncl)+1i*randn(1,Ncl))/sqrt(2);
D = gamma*diag(alpha);
% AoA = 2*pi*rand(1,Ncl)-pi;
% AoD = 2*pi*rand(1,Ncl)-pi;
% AoA = (pi/3)*rand(1,Ncl)-pi/6; % uniform in [-pi/6,pi/6]
% AoD = (pi/3)*rand(1,Ncl)-pi/6; % uniform in [-pi/6,pi/6]
% AoA = pi*rand(1,Ncl)-pi/2; % uniform in [-pi/2,pi/2]
% AoD = pi*rand(1,Ncl)-pi/2;% uniform in [-pi/2,pi/2]
AoA = pi*rand(1,Ncl); % uniform in [0,2*pi]
AoD = pi*rand(1,Ncl); % uniform in [0,2*pi]
d = lamada/2; 
for l=1:Ncl
    ABS(:,l) = array_respones(AoD(l),Nt,d,lamada);
    AMS(:,l) = array_respones(AoA(l),Nr,d,lamada);
end
H = AMS*D*ABS';
[U,~,V] = svd(H);
Fopt = V(:,1:Ns);
d = diag(D);
Wopt = U(:,1:Ns);
% Wopt = (1/sqrt(SNR) * inv(Fopt'*(H'*H)*Fopt + 1/SNR*eye(Ns))*Fopt'*H')'; % rx, mmse
end