function [FRF,S] = FastAntSelTX(ABS,pu,H,NRF)
Nt = size(H,2); %#ok<*NASGU>
Nr = size(H,1);
Nps=size(ABS,2);
S = [];
I = 1:Nps;
G = eye(Nr);
g = zeros(1,Nps);
for j = 1:Nps
    g(j) = ABS(:,j)'*(H'*H)*ABS(:,j);
end
for n = 1:NRF
    [~,J] = max(abs(g));
    idx = (I == J); 
    I(idx) = [];
    S = [S,J];
    if n < NRF
        a_vec = 1/sqrt(1/pu + g(J))*G*H*ABS(:,J);
        G = G - a_vec*a_vec';
        g(J) = 0;
        for j = 1:Nps-n
            g(I(j)) = g(I(j)) - abs(ABS(:,I(j))'*H'*a_vec)^2;
        end
    end
end
FRF = ABS(:,S);
