function [WRF,S] = FastAntSelRX(AMS,pu,He,NRFt, NRFr)
Nr = size(He,1); %#ok<*NASGU>
Nps=size(AMS,2);
S = [];
I = 1:Nps;
G = eye(NRFt);
g = zeros(1,Nps);
for j = 1:Nps
    g(j) = AMS(:,j)'*(He*He')*AMS(:,j);
end
for n = 1:NRFr
    [~,J] = max(abs(g));
    idx = (I == J); 
    I(idx) = [];
    S = [S,J];
    if n < NRFr
        a_vec = 1/sqrt(1/pu + g(J))*G*He'*AMS(:,J);
        G = G - a_vec*a_vec';
        g(J) = 0;
        for j = 1:Nps-n
            g(I(j)) = g(I(j)) - abs(AMS(:,I(j))'*He*a_vec)^2;
        end
    end
end
WRF = AMS(:,S);
