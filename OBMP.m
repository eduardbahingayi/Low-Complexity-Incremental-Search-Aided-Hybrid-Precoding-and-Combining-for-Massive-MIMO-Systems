function [ FRF, FBB ] = OBMP(Fopt, At, NRF, Ns )
%Orthogonal Based Matching Pursuit
PU = At' * Fopt;
B = diag(PU*PU');
S = [];
for n = 1:NRF
    [~,k] = max(B);
    S = [S,k];
    B(k) = 0;
end
    FRF = At(:,S);
    FBB = PU(S,:);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
end

