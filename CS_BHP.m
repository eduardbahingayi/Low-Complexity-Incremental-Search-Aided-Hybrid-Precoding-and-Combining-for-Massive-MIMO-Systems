function [FRF, FBB] = CS_BHP(Hb_bar,ADFT,NRF,Ns)
% Function to generate Vb;
[~,Vb] = BaemSpaceSVD(Hb_bar,ADFT,Ns);
gamma = diag(Vb*Vb');
W = [];
for l = 1:NRF
    [~,k] = max(gamma);
    W = [W,k];
    gamma(k) = 0;
end
FRF = ADFT(:,W);
FBB = Vb(W,1:Ns);
FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');

