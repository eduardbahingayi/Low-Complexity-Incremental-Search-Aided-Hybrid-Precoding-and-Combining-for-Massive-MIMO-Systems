function [AtDFT] = DFT_Codebook(N)
Ant_index = (1:N)-1;
Dict_Res  = (0:1:N-1)/N;
AtDFT     = sqrt(1/N)*exp(1i*(2*pi)*Ant_index'*(Dict_Res));