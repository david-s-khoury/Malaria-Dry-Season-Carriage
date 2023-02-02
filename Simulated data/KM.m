function [t,n,m,q,S] = KM(times, t_max)
% Kaplan-Meier curves
% Input:
%   times   vector with survival times (without 0, only censored if times>t_max)
%   t_max   maximal time (i.e., surveillance period)
% Output:
%   t       vector with ordered and unique survival times
%   n       number of subjects in the risk set (same length as t)
%   m       number of events  
%   q       number of censored
%   S       survival probability

tmp = [0,sort(times(times<=t_max))];
t = unique(tmp);
n = [sum(times>0);zeros(length(t)-1,1)]; 
m = zeros(length(t),1);
q = zeros(length(t),1);
S = [1;zeros(length(t)-1,1)];

for i = 2:length(t)
    if i==length(t) && any(times>t_max)
        q(i) = sum(times>t_max);
    end
    n(i) = n(i-1) - m(i-1) - q(i-1);
    m(i) = sum(times==t(i));
    S(i) = S(i-1)*(n(i) - m(i))/n(i);
end

if any(times>t_max)
    t = [t,t_max];
    S = [S;S(length(S))];
end

end

