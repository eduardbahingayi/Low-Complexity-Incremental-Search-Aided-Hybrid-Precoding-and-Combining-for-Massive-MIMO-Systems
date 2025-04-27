function [ FRF, FBB ] = HybridlyConnectedTX(H,N,Ns)
Q = H'*H;
[U,E,V] = svd(Q);
F = U(:,1:Ns);
B = (F'*F)^(1/2);
P = F*inv(B);
FRF = sqrt(1/N)*exp (1i*angle(P));
FBB = pinv(FRF)*F;
   