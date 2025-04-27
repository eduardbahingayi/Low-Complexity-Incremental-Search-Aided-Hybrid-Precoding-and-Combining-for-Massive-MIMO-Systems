function [ FRF, FBB ] = OMP(Fopt, At, NRF )
FRF = [];
Fres = Fopt;
S = [];
for n = 1:NRF
    PU = At' * Fres;
    [~,k] = max(sum( abs(PU).^2, 2 ));
    S=[S,k];
    FRF = [FRF , At(:,k)];
    FBB = pinv(FRF) * Fopt; %use pseudoinverse to avoid the inverse of a possible singular matrix
    Fres = (Fopt - FRF * FBB) / norm(Fopt - FRF * FBB,'fro');
end
end

