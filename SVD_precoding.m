function [F,W] = SVD_precoding(H,Ns)
[U,~,V]=svd(H);
F = V(:,1:Ns);
W = U(:,1:Ns);



