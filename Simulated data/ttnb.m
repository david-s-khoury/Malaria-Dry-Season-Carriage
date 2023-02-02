function [time] = ttnb(ct,bites_dry,bites_mal)
% Computes the time to the next mosquito bite
% Input:
%   ct          current time (in the year)
%   bites_dry   biting rate during the dry season
%   bites_mal   biting rate during the malaria transmission season
% Output:
%   time        time to next bite

% varying biting rate depending on the season:
bite = @(t) interp1([1,181,182,365,366],[bites_dry,bites_dry,bites_mal,bites_mal,bites_dry],t);

day = mod(floor(ct)-1,365)+1; % day of the year

% We use a non-homogeneous Poisson-process to compute the time to the next
% bite (see Supplementary methods for details):
x = -log(1-rand(1));
tmp = (ceil(ct)-ct)*bite(day);
count = 0;
while(tmp<x)
    count = count+1;
    tmp = tmp+bite(mod(day+count-1,365)+1);
end
d = (tmp-x)/bite(mod(day+count-1,365)+1);

if count==0
    time = (ceil(ct)-ct) - d;
else
    time = (ceil(ct)-ct) + count-1 + (1-d);
end

end