function[FBB,WBB] = Baseband(WRF,H,FRF,Ns) 
%This function design digital precoders/combiners
%Ns  = Number of Data streams
%NRF = Number of RF chains
%FRF = Analog Precoder
%WRF = Analog combiner
%H   = Channel matrix
Heq     = H*FRF;
[U,~,V] = svd(WRF'*Heq);
FBB     = (FRF'*FRF)^(-0.5)*V(:,1:Ns);
FBB     = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
WBB     = U(:,1:Ns);       
end
  