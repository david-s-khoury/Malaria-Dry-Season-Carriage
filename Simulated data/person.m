function [P,S,C] = person(age,birth,bites_dry,bites_mal,par)
% simulates the parasite concentration P, strain specific immunity S, and
% cross-reactive immunity C of a person with a certain age
% includes biting rate varying by season and for each individual, time to
% next bite computed using ttnb.m, and each bite inoculates a new strain
% intra-host dynamics using intrahost.m

% Input:
%   age         age of the person (in days)
%   birth       day of birth (in the year)
%   bites_dry   biting rate during the dry season
%   bites_mal   biting rate during the malaria transmission season
%   par         parameters of the model:
%       Z_p     parsite clearance threshold
%       r       initial PMR 
%       a       strain specific immunity acquisition rate
%       g       cross-reactive immunity acquisition rate
%       d       cross-reactive immunity loss rate
%       e       initial number of infected RBCs (in millions)
%       liver   duration of the liver stage (in days)
% Output:
%   P           parasite concentration, size: nx1
%   S           specific immunity strength, size: nx1
%   C           cross-reactive immunity strength, size: 1x1


% blood volume (in litres) depending on age (in days)
V =@(a) (a<22*365).*(0.00059.*a+0.3) + (a>=22*365).*5;

n = 500; % number of rows in P and S

% 1st bite
tcb = ttnb(birth,bites_dry,bites_mal); % time to current (first) bite
tnb = tcb + ttnb(birth+tcb,bites_dry,bites_mal); % time of next bite

y0 = zeros(2*n+1,1);
y0(1,1) = par{1,'e'}/V(tcb);  % transfer of parasites by bite
if (tnb+par{1,'liver'})>age
    [~,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+par{1,'liver'}) age],y0);
else
    [~,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),...
        [(tcb+par{1,'liver'}) (tnb+par{1,'liver'})],y0);
end

tcb = tnb;
P = (y(:,1:n))';
S = (y(:,(n+1):(2*n)))';
C = (y(:,2*n+1))';


% later bites
while((tcb+par{1,'liver'})<age)
    tnb = tcb + ttnb(birth+tcb,bites_dry,bites_mal); % time of next bite
    
    % if parasites are under a threshold, set parasite concentration to 0
    % also set strain specific immunity to zero(it would only decay 
    % exponentially and not influence C or other parasite strains)
    if any(P(1:(find(max(P,[],2)==0,1)-1),:)<par{1,'Z_p'},'all')
        for i=1:(find(max(P,[],2)==0,1)-1)
            tmp = find(P(i,:)<par{1,'Z_p'},1);
            P(i,tmp:size(P,2)) = zeros(size(tmp:size(P,2)));
            S(i,tmp:size(P,2)) = zeros(size(tmp:size(P,2)));
        end
    end
    
    y0 = [P(:,size(P,2));S(:,size(P,2));C(size(P,2))]; % initial condition
    if any(P(:,size(P,2))==0) % find next free row for new parasite strain
        cs = find(P(:,size(P,2))==0,1);
    else
        disp('error: not enough rows in P and S')
        break
    end
    y0(cs,1) = y0(cs,1) + par{1,'e'}/V(tcb);  % transfer of parasites by bite
    
    if (tnb+par{1,'liver'})>age
        [~,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+par{1,'liver'}) age],y0);
    else
        [~,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),...
            [(tcb+par{1,'liver'}) (tnb+par{1,'liver'})],y0);
    end

    tcb = tnb;
    P = (y(:,1:n))';
    S = (y(:,(n+1):(2*n)))';
    C = (y(:,2*n+1))';
end


% if parasites are under a threshold, set parasite concentration to 0
% also set strain specific immunity to zero (it would only decay 
% exponentially and not influence C or other parasite strains)
if any(P(1:(find(max(P,[],2)==0,1)-1),:)<par{1,'Z_p'},'all')
    for i=1:(find(max(P,[],2)==0,1)-1)
        tmp = find(P(i,:)<par{1,'Z_p'},1);
        P(i,tmp:size(P,2)) = zeros(size(tmp:size(P,2)));
        S(i,tmp:size(P,2)) = zeros(size(tmp:size(P,2)));
    end
end

ind = size(P,2);
P = P(:,ind);
S = S(:,ind);
C = C(ind);

end

