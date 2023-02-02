function [out] = intrahost(P,S,C,par)
% deterministic intrahost dynamics, i.e., the right-hand side of the ODE
% Input:
%   P       parasite concentration, size: nx1 (n is the number of strains)
%   S       strain specific immunity strength, size: nx1
%   C       cross-reactive immunity strength, size: 1x1
%   par     parameter vector containing:
%       Z_p     parsite clearance threshold
%       r       initial PMR 
%       a       strain specific immunity acquisition rate
%       g       cross-reactive immunity acquisition rate
%       d       cross-reactive immunity loss rate
% Output:
%   out     right-hand side of the ODE for intrahost dynamics

n = length(P);
out_P = (log(par{1,'r'})/2-C*ones(n,1)-S).*P.*(P>=par{1,'Z_p'});
out_S = par{1,'a'}.*P;
out_C = par{1,'g'}*sum(P(P>=par{1,'Z_p'}))-par{1,'d'}*C;

out = [out_P;out_S;out_C];
