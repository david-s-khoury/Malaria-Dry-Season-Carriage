% Model simulations and plots for the article "Evidence for exposure 
% dependent carriage of malaria parasites across the dry season: modelling
% analysis of longitudinal data"

%%% This file contains the code for the following model simulations:
% 1. 1 individual from age 0 with seasonal mosquito biting rate
% 2. 1000 individuals with random mosquito biting rates
% 3. N individuals with specific biting rates


%%% This code uses the following functions:
% ttnb          to compute the time to the next bite
% intrahost     deterministic intrahost dynamics
% person        to compute the parasite density and immunity of an
%               individual of a certain age
% KM            for Kaplan-Meier curves


%%%%%%%%%%%%%%%%%%
%%% Parameters %%%
%%%%%%%%%%%%%%%%%%

r = 14; % initial PMR (per cycle)
a = 7e-6; % strain specific immunity acquisition rate
g = 1e-6; % general immunity acquisition rate
d = 3.7e-4; % general immunity loss rate
e = 0.056; % initially infected RBCs (in millions)
lod = 1; % detection threshold (parasites/microlitre)
lod_rdt = 100; % limit of detection for RDT (parasites/microliter)
Z_p = 0.005; % parasite clearance threshold (parasites/microlitre)
liver = 7; % duration of liver-stage (in days)
eds = 181; % end of dry season (day of the year), June 30th:  181st day
bites_dry = 0; % biting rate in the dry season (in bites per day)
bites_mal_min = 0.04; % minimal biting rate in the malaria transmission season (bites per day)
bites_mal_max = 0.04; % maximal biting rate in the malaria transmission season (bites per day)

% table with parameters:
par = table('Size',[1,13],'VariableTypes',{'double','double','double','double','double','double'...
    'double','double','double','double','double','double','double'},...
    'VariableNames',{'g','r','d','a','e','lod','lod_rdt','Z_p','liver','eds',...
    'bites_dry','bites_mal_min','bites_mal_max'});
par{1,'g'} = g; par{1,'r'} = r; par{1,'d'} = d; par{1,'a'} = a;
par{1,'e'} = e; par{1,'lod'} = lod; par{1,'lod_rdt'} = lod_rdt; par{1,'Z_p'} = Z_p;
par{1,'liver'} = liver; par{1,'eds'} = eds; par{1,'bites_dry'} = bites_dry;
par{1,'bites_mal_min'} = bites_mal_min; par{1,'bites_mal_max'} = bites_mal_max;

% blood volume (in litres) depending on age (in days)
V =@(a) (a<22*365).*(0.00059.*a+0.3) + (a>=22*365).*5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 1. 1 individual from age 0 with seasonal mosquito biting rate %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parasite concentration, specific and general immunity strength
n = 500; % number of strains
t_max = 365*20; % maximal simulated time interval (in days)
bt = zeros(1,100*t_max/365); % bite times
t = zeros(3100*t_max/365,1); % time in days
P = zeros(n,3100*t_max/365); % in parasites/microlitre
S = zeros(n,3100*t_max/365); % strain specific immunity
C = zeros(1,3100*t_max/365); % general immunity

%%%%%%%%%%%%%%%%
%%% 1st bite %%%
%%%%%%%%%%%%%%%%

birth = rand(1)*365+1; % random time of birth
if birth==366
    birth = 1;
end
% biting rate during the malaria transmission season:
% random rate between minimal and maximal bites per day
bites_mal_ind = bites_mal_min + rand(1)*(bites_mal_max-bites_mal_min);

tcb = ttnb(birth,bites_dry,bites_mal_ind); % time to first bite
tnb = tcb + ttnb(birth+tcb,bites_dry,bites_mal_ind); % time of next bite
bt(1) = tcb; % bite times

% deterministic intrahost dynamics:
y0 = zeros(2*n+1,1); % initial condition: no parasites, no immunity
y0(1,1) = e/V(tcb);  % transfer of parasites by bite
[tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+liver) (tnb+liver)],y0);
t(2:length(tmp),1) = tmp(1:(length(tmp)-1));
P(:,2:length(tmp)) = (y(1:(length(tmp)-1),1:n))';
S(:,2:length(tmp)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
C(1,2:length(tmp)) = (y(1:(length(tmp)-1),2*n+1))';
tcb = tnb;
ind = 1; % number of bite
ind2 = length(tmp)+1; % next free column in t, P, S and C


%%%%%%%%%%%%%%%%%%%
%%% later bites %%%
%%%%%%%%%%%%%%%%%%%

while(max(t)<t_max)
    ind = ind+1;
    bt(ind) = tcb; % save time of current bite
    tnb = tcb + ttnb(birth+tcb,bites_dry,bites_mal_ind); % time of next bite
    
    % if parasites are under a threshold, set parasite concentration to 0
    % also set strain specific immunity to zero
    % (it would only decay exponentially and not influence G or other parasite strains)
    if any(P(1:(find(max(P,[],2)==0,1)-1),find(t==bt(ind-1)+liver):(ind2-1))<Z_p,'all')
        for i=1:(find(max(P,[],2)==0,1)-1)
            tmp = find(P(i,find(t==bt(ind-1)+liver):(ind2-1))<Z_p,1);
            P(i,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
            S(i,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
        end
    end
    
    y0 = [P(:,ind2-1);S(:,ind2-1);C(1,ind2-1)]; % initial condition
    
    % put the new parasite strain in the first empty row
    if any(P(:,ind2-1)==0)
        cs = find(P(:,ind2-1)==0,1);
    else
        break
    end
    y0(cs,1) = y0(cs,1) + e/V(tcb);
    
    [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+liver) (tnb+liver)],y0);
    t(ind2:(ind2+length(tmp)-2),1) = tmp(1:(length(tmp)-1));
    P(:,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),1:n))';
    S(:,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
    C(1,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),2*n+1))';
    tcb = tnb;
    ind2 = ind2+length(tmp)-1;
end

rows = find(max(P,[],2)==0,1)-1;
bt = bt(1:ind);
t = t(1:ind2-1,1);
P = P(1:rows,1:(ind2-1));
S = S(1:rows,1:(ind2-1));
C = C(1,1:(ind2-1));

% As before, clear parasites and strain specific immunity if the parasite
% concentration is below the clearance threshold:
if any(P(:,find(t==bt(length(bt))+liver):length(t))<Z_p,'all')
    for i=1:rows
        tmp = find(P(i,find(t==bt(length(bt))+liver):length(t))<Z_p,1);
        P(i,(find(t==bt(length(bt))+liver)+tmp-1):length(t)) = ...
            zeros(size((find(t==bt(length(bt))+liver)+tmp-1):length(t)));
        S(i,(find(t==bt(length(bt))+liver)+tmp-1):length(t)) = ...
            zeros(size((find(t==bt(length(bt))+liver)+tmp-1):length(t)));
    end
end

% Duration and maximal parasite concentration of an infection for each
% infection:
duration = zeros(size(bt)); % duration of infections from the beginning of the blood-stage to clearance
max_par = zeros(size(bt));
for i=1:length(bt)
    ind1 = find(t==(bt(i)+liver));
    cs = find(P(:,ind1-1)==0,1);
    if any(P(cs,ind1:length(t))==0)
        ind2 = find(P(cs,ind1:length(t))==0,1)+ind1-1;
        duration(i) = t(ind2-1)-t(ind1);
    else
        ind2 = length(t);
        duration(i) = nan;
    end
    if any(P(cs,ind1:ind2)>lod)
        [max_tmp,ind_max] = max(P(cs,ind1:ind2));
        ind_max = ind_max+ind1-1;
        if any(P(cs,ind_max:ind2)<lod)
            max_par(i) = max_tmp;
        else
            max_par(i) = nan;
        end
    else
        [max_tmp,ind_max] = max(P(cs,ind1:ind2));
        max_par(i) = max_tmp;
    end
end


%%%%%%%%%%%%%
%%% Plots %%%
%%%%%%%%%%%%%

load('Data-1Ind-low-FOI.mat') % for Fig. S1
% load('Data-1Ind-med-FOI.mat') % for Fig. 2
% load('Data-1Ind-high-FOI.mat') % for Fig. S2

% Overall parasite concentration (log scale):
subplot(3,2,1:2)
paras = sum(P,1);
semilogy(t/365,paras,'Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0 t(length(t))/365])
ylim([5e-3,1.5e5])
hold on
xlabel('Age [years]','Fontsize',15)
ylabel('Parasites/microliter','Fontsize',15)
title('Concentration of parasites (all strains)','Fontsize',15)
hold off
%     yline(lod_rdt,'blue','LineWidth',2)
%     yline(lod,'red','LineWidth',2)

% Strain specific and cross-reactive immunity:
subplot(3,2,3:4)
plot(t/365,S,'LineWidth',2,'color',[0.5 0.5 0.5])
hold on
plot(t/365,C,'r','LineWidth',5)
yline(0,'k','LineWidth',2)
xlim([0 t(length(t))/365])
ylim([-0.2/7*3 3])
% indicate bites
for i=1:length(bt)
    triang = polyshape([-0.03 0.03 0 -0.03]+bt(i)/365,[-1,-1,0,-1]*0.2*3/7);
    h=plot(triang,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',1);
    h.EdgeColor = [0.5 0.5 0.5];
end
xlabel('Age [years]','Fontsize',15)
ylabel('Immunity level','Fontsize',15)
title('Strain specific and cross-reactive immunity','Fontsize',15)
hold off

% Duration of infections:
subplot(3,2,5)
plot(bt/365,duration,'k.','MarkerSize',15)
xlabel('Age [years]','Fontsize',15)
ylabel('Duration of an infection [days]','Fontsize',15)
title('Duration of infections','Fontsize',15)

% Peak parasite concentration:
subplot(3,2,6)
semilogy(bt/365,max_par,'k.','MarkerSize',15)
xlabel('Age [years]','Fontsize',15)
ylabel('Parasites/microliter','Fontsize',15)
title('Maximal parasite conc. during an infection','Fontsize',15)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 2. 1000 individuals with random mosquito biting rates %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate 1000 individuals with different ages at the end of the dry season
% RDT positive individuals are treated.
% In the code, PCR positives are classified as "carriers", otherwise
% they are non-carriers.
% Simulate all individuals for the duration of the surveillance period.

% for efficacy we simulate 125 times 8 individuals in parallel:
for l=1:125 % 125 times 8 individuals
    
    r = 14; % initial PMR (per cycle)
    a = 7e-6; % strain specific immunity acquisition rate
    g = 1.5e-6; % general immunity acquisition rate
    d = 3.7e-4; % general immunity loss rate
    e = 0.056; % initially infected RBCs (in millions)
    lod = 1; % detection threshold for PCR (parasites/microlitre)
    lod_rdt = 100; % limit of detection for RDT (parasites/microliter)
    Z_p = 2e-7; % parasite clearance threshold (parasites/microlitre)
    liver = 7; % duration of liver-stage (in days)
    eds = 181; % end of dry season (day of the year), June 30th:  181st day
    bites_dry = 0; % biting rate in the dry season (in bites per day)
    %%% for simulating a heterogeneous popualtions with different FOIs:
    %     bites_mal_min = 1/250; % minimal biting rate in the malaria transmission season (bites per day)
    %     bites_mal_max = 1/25; % maximal biting rate in the malaria transmission season (bites per day)
    %%% for simulating a homogeneous population with the same FOI for all:
    bites_mal_min = 0.04; % minimal biting rate in the malaria transmission season (bites per day)
    bites_mal_max = 0.04; % maximal biting rate in the malaria transmission season (bites per day)
    
    % table with parameters:
    par = table('Size',[1,13],'VariableTypes',{'double','double','double','double','double'...
        'double','double','double','double','double','double','double','double'},...
        'VariableNames',{'g','r','d','a','e','lod','lod_rdt','Z_p','liver','eds',...
        'bites_dry','bites_mal_min','bites_mal_max'});
    par{1,'g'} = g; par{1,'r'} = r; par{1,'d'} = d; par{1,'a'} = a;
    par{1,'e'} = e; par{1,'lod'} = lod; par{1,'lod_rdt'} = lod_rdt; par{1,'Z_p'} = Z_p;
    par{1,'liver'} = liver; par{1,'eds'} = eds; par{1,'bites_dry'} = bites_dry;
    par{1,'bites_mal_min'} = bites_mal_min; par{1,'bites_mal_max'} = bites_mal_max;
    
    % blood volume (in litres) depending on age (in days)
    V =@(a) (a<22*365).*(0.00059.*a+0.3) + (a>=22*365).*5;
    
    
    N = 8; % number of individuals
    n = 500; % number of strain (rows)
    age_min = 0.25*365; % minimal age at the end of the dry season (days)
    age_max = 12*365; % maximal age at the end of the dry season (days)
    surv = 1*365; % surveillance period (days)
    
    carrier = zeros(1,N); % 1 if individual is carrier
    treated = zeros(1,N); % 1 if individual is treated
    ages = zeros(N,1); % ages at the beginning of treatment
    A_tmp = cell(N,1); % save age, PMR, t, P, S, C, bt, and biting rate for each individual
    
    parfor i=1:N
        age = age_min + rand(1)*(age_max-age_min); % random age at the end of the dry season
        birth = mod(eds-1-age,365)+1; % day of birth (in the year)
        bites_mal_ind = bites_mal_min + rand(1)*(bites_mal_max-bites_mal_min); % biting rate during the transmission season
        % simulate a person from birth to their age:
        [P_tmp,S_tmp,C_tmp] = person(age,birth,bites_dry,bites_mal_ind,par); % parasites & immunity of person
        
        bt = zeros(1,100*surv/365); % bite times
        t = zeros(1+3100*surv/365,1); % time in days
        P = [P_tmp,zeros(n,3100*surv/365)]; % parasite concentration
        S = [S_tmp,zeros(n,3100*surv/365)]; % strain specific immunity
        C = [C_tmp,zeros(1,3100*surv/365)]; % general immunity
        
        carrier(1,i) = sum(P_tmp)>lod; % carrier if PCR positive
        treated(1,i) = sum(P_tmp)>lod_rdt; % treated if RDT positive
        
        % surveillance:
        tfb = ttnb(eds,bites_dry,bites_mal_ind); % time to first bite
        if treated(1,i)==1 % case with treatment
            % time to first bite:
            y0 = [P(:,2);S(:,2);C_tmp]; % the second row of P and S consists of 0s only
            if (tfb+liver)>surv
                [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),...
                    y(2*n+1,1),par),[0 surv],y0);
                t(2:(length(tmp)+1),1) = tmp(:,1);
                P(:,2:(length(tmp)+1)) = (y(:,1:n))';
                S(:,2:(length(tmp)+1)) = (y(:,(n+1):(2*n)))';
                C(1,2:(length(tmp)+1)) = (y(:,2*n+1))';
                ind2 = length(tmp)+2;
            else
                [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),...
                    y(2*n+1,1),par),[0 (tfb+liver)],y0);
                t(2:length(tmp),1) = tmp(1:(length(tmp)-1),1);
                P(:,2:length(tmp)) = (y(1:(length(tmp)-1),1:n))';
                S(:,2:length(tmp)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
                C(1,2:length(tmp)) = (y(1:(length(tmp)-1),2*n+1))';
                ind2 = length(tmp)+1;
            end
        else % case without treatment:
            % time to first bite:
            y0 = [P_tmp(:,1);S_tmp(:,1);C_tmp(1)];
            if (tfb+liver)>surv
                [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),...
                    y(2*n+1,1),par),[0 surv],y0);
                t(1:length(tmp),1) = tmp(:,1);
                P(:,1:length(tmp)) = (y(:,1:n))';
                S(:,1:length(tmp)) = (y(:,(n+1):(2*n)))';
                C(1,1:length(tmp)) = (y(:,2*n+1))';
                ind2 = length(tmp)+1;
            else
                [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),...
                    y(2*n+1,1),par),[0 (tfb+liver)],y0);
                t(1:(length(tmp)-1),1) = tmp(1:(length(tmp)-1),1);
                P(:,1:(length(tmp)-1)) = (y(1:(length(tmp)-1),1:n))';
                S(:,1:(length(tmp)-1)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
                C(1,1:(length(tmp)-1)) = (y(1:(length(tmp)-1),2*n+1))';
                ind2 = length(tmp);
            end
        end
        tcb = tfb;
        ind = 0;
        
        % next bites:
        while((tcb+liver)<surv)
            ind = ind+1;
            bt(ind) = tcb; % save time of current bite
            tnb = tcb + ttnb(eds+tcb,bites_dry,bites_mal_ind); % time of next bite
            
            % if parasites are under a threshold, set parasite concentration to 0
            % also set strain specific immunity to zero (it would only decay 
            % exponentially and not influence G or other parasite strains)
            if any(P>0,'all')
                if ind==1
                    if any(P(max(P,[],2)>0,1:(ind2-1))<Z_p,'all')
                        k_ind = find(max(P,[],2)>0);
                        for j=1:(length(k_ind))
                            k = k_ind(j);
                            tmp = find(P(k,1:(ind2-1))<Z_p,1);
                            P(k,tmp:(ind2-1)) = zeros(size(tmp:(ind2-1)));
                            S(k,tmp:(ind2-1)) = zeros(size(tmp:(ind2-1)));
                        end
                    end
                else
                    if any(P(max(P,[],2)>0,find(t==bt(ind-1)+liver):(ind2-1))<Z_p,'all')
                        k_ind = find(max(P,[],2)>0);
                        for j=1:(length(k_ind))
                            k = k_ind(j);
                            tmp = find(P(k,find(t==bt(ind-1)+liver):(ind2-1))<Z_p,1);
                            P(k,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = ...
                                zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
                            S(k,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = ...
                                zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
                        end
                    end
                end
            end
            
            y0 = [P(:,ind2-1);S(:,ind2-1);C(1,ind2-1)];
            % save parasite concentration in first "empty" row:
            cs = 1;
            if any(P(:,ind2-1)==0)
                cs = find(P(:,ind2-1)==0,1);
            else
                disp('error: not enough rows')
                break
            end
            y0(cs,1) = y0(cs,1) + e/V(tcb);
            if (tnb+liver)>surv
                [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),...
                    y(2*n+1,1),par),[(tcb+liver) surv],y0);
                t(ind2:(ind2+length(tmp)-1),1) = tmp(1:length(tmp));
                P(:,ind2:(ind2+length(tmp)-1)) = (y(1:length(tmp),1:n))';
                S(:,ind2:(ind2+length(tmp)-1)) = (y(1:length(tmp),(n+1):(2*n)))';
                C(1,ind2:(ind2+length(tmp)-1)) = (y(1:length(tmp),2*n+1))';
            else
                [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),...
                    y(2*n+1,1),par),[(tcb+liver) (tnb+liver)],y0);
                t(ind2:(ind2+length(tmp)-2),1) = tmp(1:(length(tmp)-1));
                P(:,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),1:n))';
                S(:,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
                C(1,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),2*n+1))';
            end
            tcb = tnb;
            ind2 = ind2+length(tmp)-1;
        end
        
        rows_ind = max(P,[],2)>0;
        if ind==0
            bt = [];
        else
            bt = bt(1:ind);
        end
        if any(rows_ind==1)
            t = t(1:ind2-1,1);
            P = P(rows_ind,1:(ind2-1));
            S = S(rows_ind,1:(ind2-1));
            C = C(1,1:(ind2-1));
        else
            t = t(1:ind2-1,1);
            P = P(1,1:(ind2-1));
            S = S(1,1:(ind2-1));
            C = C(1,1:(ind2-1));
        end
        
        % if parasites are under a threshold, set parasite concentration to 0
        % also set strain specific immunity to zero
        % (it would only decay exponentially and not influence G or other parasite strains)
        if ind==0
            if any(P<Z_p,'all')
                for j=1:size(P,1)
                    tmp = find(P(j,:)<Z_p,1);
                    P(j,tmp:length(t)) = zeros(size(tmp:length(t)));
                    S(j,tmp:length(t)) = zeros(size(tmp:length(t)));
                end
            end
        else
            if any(P(:,find(t==bt(length(bt))+liver):length(t))<Z_p,'all')
                for j=1:size(P,1)
                    tmp = find(P(j,find(t==bt(length(bt))+liver):length(t))<Z_p,1);
                    P(j,(find(t==bt(length(bt))+liver)+tmp-1):length(t)) = ...
                        zeros(size((find(t==bt(length(bt))+liver)+tmp-1):length(t)));
                    S(j,(find(t==bt(length(bt))+liver)+tmp-1):length(t)) = ...
                        zeros(size((find(t==bt(length(bt))+liver)+tmp-1):length(t)));
                end
            end
        end
        
        % PMR for each strain at the end of the surveillance period:
        PMR = r.*exp(-2.*(S(:,size(S,2))+C(length(C)).*ones(size(S,1),1)));
        
        ages(i) = age;
        A_tmp{i,1}={age,PMR,t,P,S,C,bt,bites_mal_ind};
        
        % progress:
        %         disp(i)
    end
    
    A = cell(N,8);
    for i=1:N
        for j=1:8
            A{i,j} = A_tmp{i}{j};
        end
    end
    
    filename = ['Data-prel-' num2str(l) '.mat'];
    save(filename,'A','par')
    
    % progress:
    % disp('.')
    
    clear all
    
end

% put together the l files into one file:
A_tmp = cell(100,8);

for l=1:125
    filename = ['Data-prel-' num2str(l) '.mat'];
    load(filename)
    
    for i=1:8
        for j=1:8
            A_tmp{i+(l-1)*8,j} = A{i,j};
        end
    end
    
end

A = A_tmp;
%%% for simulating a heterogeneous popualtions with different FOIs:
% save('Data-1','A','par')
%%% for simulating a homogeneous population with the same FOI for all:
% save('Data-3-10','A','par')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots and analysis of simulation of 1,000 individuals with heterogeneous FOI and age %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
load('Data-1.mat')

% load parameters:
% r = par{1,'r'}; d = par{1,'d'}; a = par{1,'a'}; b = par{1,'b'}; g = par{1,'g'};
% e = par{1,'e'}; lod = par{1,'lod'}; lod_rdt = par{1,'lod_rdt'}; Z_p = par{1,'Z_p'}; Z_s = par{1,'Z_s'};
% Z_g = par{1,'Z_g'}; liver = par{1,'liver'}; eds = par{1,'eds'}; bites_dry = par{1,'bites_dry'};
bites_mal_min = par{1,'bites_mal_min'}; bites_mal_max = par{1,'bites_mal_max'};
% gmax = par{1,'gmax'};

ages = horzcat(A{:,1})';
mean_bites = horzcat(A{:,8});
N = length(ages);
carrier = zeros(N,1);
treated = zeros(N,1);
for i = 1:N
    P_tmp = A{i,4};
    carrier(i,1) = sum(P_tmp(:,1))>par{1,'lod'};
    treated(i,1) = sum(P_tmp(:,1))>par{1,'lod_rdt'};
end


%%% FOI and age heatmap: (Fig. S9b)
m_bites = vertcat(A{:,8});
age_groups = [0,365,2*365,3*365,4*365,5*365,6*365,7*365,8*365,9*365,10*365,11*365,max(ages)+1];
age_groups_names = {'<1','1','2','3','4','5','6','7','8','9','10','>=11'};
age_group = strings(N,1);
% exp_groups = [0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02];
% exp_groups_names = {'0.004-0.006','0.006-0.008','0.008-0.01','0.01-0.012','0.012-0.014',...
%     '0.014-0.016','0.016-0.018','0.018-0.02'};
exp_groups = [0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032,0.036,0.04];
exp_groups_names = {'4-8','8-12','12-16','16-20','20-24','24-28','28-32','32-36','36-40'};
% exp_groups = [0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05];
% exp_groups_names = {'5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50'};
m_bites_group = strings(N,1);
for i=1:N
    age_group(i) = age_groups_names{find((age_groups-ages(i))>0,1)-1};
    m_bites_group(i) = exp_groups_names{find((exp_groups-m_bites(i))>0,1)-1};
end
PCR_pos = carrier;
RDT_pos = treated;
tbl = table(ages,age_group,m_bites,m_bites_group,PCR_pos,RDT_pos);

heatmap(tbl,'age_group','m_bites_group','ColorVariable','PCR_pos');
xlabel('Age group [years]')
ylabel('Force Of Infection (FOI) [bites/1000 days]')
title('Fraction of PCR positive by age and exposure')


%%% Cross-reactive immunity by age and FOI (heatmap): (Fig. S9a)
age_groups = [0,365,2*365,3*365,4*365,5*365,6*365,7*365,8*365,9*365,10*365,11*365,max(ages)+1];
age_groups_names = {'<1','1','2','3','4','5','6','7','8','9','10','>=11'};
age_group = strings(N,1);
PCR_pos = carrier;
RDT_pos = treated;
gen_im = zeros(N,1);
gen_im_group1 = strings(N,1);
gen_im_group2 = strings(N,1);
gen_im_groups1 = [-1,0.28,0.42,0.56,0.7,0.84,0.98,1.12,1.26,2];
gen_im_group_names1 = {'<0.28','0.28-0.42','0.42-0.56','0.56-0.7','0.7-0.84','0.84-0.98','0.98-1.12',...
    '1.12-1.26','>=1.26'};
gen_im_groups2 = [-1,0.67,0.878,1.013,1.115,1.193,1.266,1.270205,1.27053,1.271045,2];
gen_im_group_names2 = {'<0.67','0.67-0.878','0.878-1.013','1.013-1.115','1.115-1.193','1.193-1.266',...
    '1.266-1.270205','1.270205-1.27053','1.27053-1.271045','>=1.271055'};
% exp_groups = [0.004,0.006,0.008,0.01,0.012,0.014,0.016,0.018,0.02];
% exp_groups_names = {'0.004-0.006','0.006-0.008','0.008-0.01','0.01-0.012','0.012-0.014',...
%     '0.014-0.016','0.016-0.018','0.018-0.02'};
exp_groups = [0.004,0.008,0.012,0.016,0.02,0.024,0.028,0.032,0.036,0.04];
exp_groups_names = {'4-8','8-12','12-16','16-20','20-24','24-28','28-32','32-36','36-40'};
% exp_groups = [0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05];
% exp_groups_names = {'5-10','10-15','15-20','20-25','25-30','30-35','35-40','40-45','45-50'};
m_bites_group = strings(N,1);
for i=1:N
    age_group(i) = age_groups_names{find((age_groups-ages(i))>0,1)-1};
    gen_im(i) = A{i,6}(1);
    gen_im_group1(i) = gen_im_group_names1{find((gen_im_groups1-gen_im(i))>0,1)-1};
    gen_im_group2(i) = gen_im_group_names2{find((gen_im_groups1-gen_im(i))>0,1)-1};
    m_bites_group(i) = exp_groups_names{find((exp_groups-m_bites(i))>0,1)-1};
end
tbl = table(ages,age_group,gen_im,gen_im_group1,gen_im_group2,m_bites_group,PCR_pos,RDT_pos);

heatmap(tbl,'age_group','m_bites_group','ColorVariable','gen_im');%,'ColorMethod','median');
xlabel('Age [years]')
ylabel('Force Of Infection (FOI) [bites/1000 days]')
title('Mean immunity by FOI and age')
colormap(parula(30))


%%% Summary of simulation results: (Fig. S8)
subplot(3,2,1) % Time from first follow-up to PCR+
% The survival curve in the final version of the manuscript was produced in
% R, see Simulated_data_analysis.R.
day = 206; % July 25th, ~ day of first follow-up
PCR_pos = zeros(1,N);
thresh = par{1,'lod'}; % use lod of PCR as the threshold
for i=1:N
    all_par = sum(A{i,4},1);
    ind = find(A{i,3}>=day-par{1,'eds'},1);
    if any(all_par(ind:length(all_par))>thresh) && all_par(ind)<thresh
        if treated(i)==1
            if any(all_par(2:size(all_par,2))>thresh)
                PCR_pos(i) = A{i,3}(ind-1+find(all_par(ind:size(all_par,2))>thresh,1))-A{i,3}(ind);
            else
                PCR_pos(i) = 1000;
            end
        else
            PCR_pos(i) = A{i,3}(ind-1+find(all_par(ind:length(all_par))>thresh,1))-A{i,3}(ind);
        end
    else
        PCR_pos(i) = 1000;
    end
end
PCR_pos_rp = PCR_pos(logical(treated));
PCR_pos_rnpn = PCR_pos(logical(1-carrier));
PCR_pos_rp = sort(PCR_pos_rp);
PCR_pos_rnpn = sort(PCR_pos_rnpn);
[t_rp,n_rp,m_rp,~,KM_rp] = KM(PCR_pos_rp,365);
[t_rnpn,n_rnpn,m_rnpn,~,KM_rnpn] = KM(PCR_pos_rnpn,365);
stairs(t_rp,KM_rp,'c','LineWidth',2)
hold on
stairs(t_rnpn,KM_rnpn,'r','LineWidth',2)
xlim([0 365])
xlabel('Time [days]','FontSize',15)
ylabel('Fraction without infection','FontSize',15)
title('Time to first detectable infection','FontSize',17)
legend({'Carriers (treated)','Non-carriers'},'Location','northeast','Fontsize',15)
hold off
subplot(3,2,2) % Mean biting rate distribution
m_bites = vertcat(A{:,8});
m_bites_tr = m_bites(logical(treated));
m_bites_non_car = m_bites(logical(1-carrier));
boxp = boxplot([m_bites;m_bites_tr;m_bites_non_car],...
    [zeros(length(m_bites),1);ones(length(m_bites_tr),1);2*ones(length(m_bites_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Mean biting rate [per day]','FontSize',15)
title('Mean biting rate distribution','FontSize',17)
[~,p_bite,ci_bite,stats_bite] = ttest2(m_bites_tr,m_bites_non_car,'tail','right');
[p_bite2,~,stats_bite2] = ranksum(m_bites_tr,m_bites_non_car,'tail','right');
subplot(3,2,3) % Fraction of RDT+ and PCR+ by age
age_groups = [0,365,2*365,3*365,4*365,5*365,6*365,7*365,8*365,9*365,10*365,11*365,max(ages)+1];
tmp_all = zeros(length(age_groups)-1,1);
tmp_pcr = zeros(length(age_groups)-1,1);
tmp_rdt = zeros(length(age_groups)-1,1);
for i=1:N
    ind = find((age_groups-ages(i))>0,1)-1;
    tmp_all(ind) = tmp_all(ind)+1;
    tmp_pcr(ind) = tmp_pcr(ind)+carrier(i)-treated(i);
    tmp_rdt(ind) = tmp_rdt(ind)+treated(i);
    if ind==1 && carrier(i)==1
        disp(i)
    end
end
ages_PCR_pos = tmp_pcr./tmp_all;
ages_RDT_pos = tmp_rdt./tmp_all;
age_groups_names = {'<1','1','2','3','4','5','6','7','8','9','10','>=11'};
bar([ages_RDT_pos';(ages_PCR_pos+ages_RDT_pos)']')
ylim([0,0.7])
ylabel('Fraction','FontSize',15)
title('Parasite carriage by age','FontSize',18)
set(gca,'XtickLabel',age_groups_names);
xlabel('Age [years]','FontSize',15);
legend({'positive RDT','positive PCR'},'Location','northwest','Fontsize',15)
subplot(3,2,4) % Age distribution of carriers and non-carriers
boxp = boxplot([ages/365;ages(logical(treated))/365;ages(logical(1-carrier))/365],...
    [zeros(length(ages),1);ones(sum(treated),1);2*ones(sum(1-carrier),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Age [years]','FontSize',15)
title('Age distribution','FontSize',17)
[p_ages,~,stats_ages] = ranksum(ages(logical(treated))/365,ages(logical(1-carrier))/365,'tail','right');
subplot(3,2,5) % PMRs for carriers and non-carriers
PMRs_all = vertcat(A{:,2})';
PMRs_all = PMRs_all(:);
PMRs_tr = vertcat(A{logical(treated),2});
PMRs_non_car = vertcat(A{logical(1-carrier),2});
boxp = boxplot([PMRs_all;PMRs_tr;PMRs_non_car],...
    [zeros(length(PMRs_all),1);ones(length(PMRs_tr),1);2*ones(length(PMRs_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('PMR','FontSize',15)
title('Parasite Multiplication Rate distribution','FontSize',17)
[~,p_pmr,ci_pmr,stats_pmr] = ttest2(PMRs_tr,PMRs_non_car,'tail','left');
[p_pmr2,~,stats_pmr2] = ranksum(PMRs_tr,PMRs_non_car,'tail','left');
subplot(3,2,6) % General immunity for carriers and non-carriers
imm_all = zeros(N,1);
for i=1:N
    imm_all(i) = A{i,6}(length(A{i,6}));
end
imm_car = imm_all(logical(treated));
imm_non_car = imm_all(logical(1-carrier));
boxp = boxplot([imm_all;imm_car;imm_non_car],...
    [zeros(length(imm_all),1);ones(length(imm_car),1);2*ones(length(imm_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Immunity','FontSize',15)
title('Immunity distribution','FontSize',17)
[p_imm,~,stats_imm] = ranksum(imm_car,imm_non_car,'tail','right');


%%% Change of infection status after 1 year of surveillance:
carrier2 = zeros(length(ages),1);
treated2 = zeros(length(ages),1);
for i=1:length(ages)
    carrier2(i,1) = sum(A{i,4}(:,size( A{i,4},2)))>par{1,'lod'};
    treated2(i,1) = sum(A{i,4}(:,size( A{i,4},2)))>par{1,'lod_rdt'};
end
% number of PCR- initially:
sum(1-carrier)
% number of RDT-PCR+ initially
sum(treated==0 & carrier==1)
% number of RDT+ initially:
sum(treated)

% transitions:
% PCR- -> PCR-:
sum(carrier==0 & carrier2==0)
% PCR- -> RDT-PCR+:
sum(carrier==0 & carrier2==1 & treated2==0)
% PCR- -> RDT+:
sum(carrier==0 & treated2==1)
% RDT-PCR+ -> PCR-:
sum(carrier==1 & treated==0 & carrier2==0)
% RDT-PCR+ -> RDT-PCR+:
sum(carrier==1 & treated==0 & carrier2==1 & treated2==0)
% RDT-PCR+ -> RDT+:
sum(carrier==1 & treated==0 & treated2==1)
% RDT+ -> PCR-:
sum(treated==1 & carrier2==0)
% RDT+ -> RDT-PCR+:
sum(treated==1 & carrier2==1 & treated2==0)
% RDT+ -> RDT+:
sum(treated==1 & treated2==1)

% number of PCR- after 1 year:
sum(1-carrier2)
% number of RDT-PCR+ after 1 year:
sum(treated2==0 & carrier2==1)
% number of RDT+ after 1 year:
sum(treated2)

% Make a Sankey diagram using R (Fig. S10, see "Simulated_data_analysis.R")


%%% Summary of simulation results where carriers/treated are those who are
%%% PCR/RDT positive after the 1 year surveillance period: (Fig. S11)
subplot(2,2,1)
age_groups = [0,365,2*365,3*365,4*365,5*365,6*365,7*365,8*365,9*365,10*365,11*365,max(ages)+1];
tmp_all = zeros(length(age_groups)-1,1);
tmp_pcr = zeros(length(age_groups)-1,1);
tmp_rdt = zeros(length(age_groups)-1,1);
for i=1:N
    ind = find((age_groups-ages(i))>0,1)-1;
    tmp_all(ind) = tmp_all(ind)+1;
    tmp_pcr(ind) = tmp_pcr(ind)+carrier2(i)-treated2(i);
    tmp_rdt(ind) = tmp_rdt(ind)+treated2(i);
    if ind==1 && carrier2(i)==1
        disp(i)
    end
end
ages_PCR_pos2 = tmp_pcr./tmp_all;
ages_RDT_pos2 = tmp_rdt./tmp_all;
age_groups_names = {'<2','2','3','4','5','6','7','8','9','10','11','>=12'};
bar([ages_RDT_pos2';(ages_PCR_pos2+ages_RDT_pos2)']')
ylim([0,0.7])
ylabel('Fraction','FontSize',15)
title('Fraction of RDT+ and PCR+ by age','FontSize',18)
set(gca,'XtickLabel',age_groups_names);
xlabel('Age [years]','FontSize',15);
legend({'positive RDT','positive PCR'},'Location','northwest','Fontsize',15)
subplot(2,2,2)
% children have aged 1 year, thus add 1 year to their age
boxp = boxplot([ages/365+1;ages(logical(treated2))/365+1;ages(logical(1-carrier2))/365+1],...
    [zeros(length(ages),1);ones(sum(treated2),1);2*ones(sum(1-carrier2),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Age [years]','FontSize',15)
title('Age distribution','FontSize',17)
[~,p_ages_2,ci_ages_2,stats_ages_2] = ttest2(ages(logical(treated2))/365,ages(logical(1-carrier2))/365,'tail',...
    'right');
[p_ages2_2,~,stats_ages2_2] = ranksum(ages(logical(treated2))/365,ages(logical(1-carrier2))/365,'tail','right');
subplot(2,2,3)
m_bites = vertcat(A{:,8});
m_bites_car = m_bites(logical(treated2));
m_bites_non_car = m_bites(logical(1-carrier2));
boxp = boxplot([m_bites;m_bites_car;m_bites_non_car],...
    [zeros(length(m_bites),1);ones(length(m_bites_car),1);2*ones(length(m_bites_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Mean biting rate [per day]','FontSize',15)
title('Mean biting rate distribution','FontSize',17)
[~,p_bite_2,ci_bite_2,stats_bite_2] = ttest2(m_bites_car,m_bites_non_car,'tail','right');
[p_bite2_2,~,stats_bite2_2] = ranksum(m_bites_car,m_bites_non_car,'tail','right');
subplot(2,2,4)
PMRs_all = vertcat(A{:,2})';
PMRs_all = PMRs_all(:);
PMRs_car = vertcat(A{logical(treated2),2});
PMRs_non_car = vertcat(A{logical(1-carrier2),2});
boxp = boxplot([PMRs_all;PMRs_car;PMRs_non_car],...
    [zeros(length(PMRs_all),1);ones(length(PMRs_car),1);2*ones(length(PMRs_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('PMR','FontSize',15)
title('Parasite Multiplication Rate distribution','FontSize',17)
[~,p_pmr_2,ci_pmr_2,stats_pmr_2] = ttest2(PMRs_car,PMRs_non_car,'tail','left');
[p_pmr2_2,~,stats_pmr2_2] = ranksum(PMRs_car,PMRs_non_car,'tail','left');


%%% Comparison of homogeneous and heterogeneous infection risk,
%%% visualizations for heterogeneous infection risk (Fig. 5d, f and h)

% Age distribution: (Fig. 5d)
boxp = boxplot([ages/365;ages(logical(treated))/365;ages(logical(1-carrier))/365],...
    [zeros(length(ages),1);ones(sum(treated),1);2*ones(sum(1-carrier),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Age [years]','FontSize',15)
title('Age distribution','FontSize',17)
[p_ages,~,stats_ages] = ranksum(ages(logical(treated))/365,ages(logical(1-carrier))/365,'tail','right');

% Immunity distribution: (Fig. 5f)
imm_all = zeros(N,1);
for i=1:N
    imm_all(i) = A{i,6}(length(A{i,6}));
end
imm_car = imm_all(logical(treated));
imm_non_car = imm_all(logical(1-carrier));
boxp = boxplot([imm_all;imm_car;imm_non_car],...
    [zeros(length(imm_all),1);ones(length(imm_car),1);2*ones(length(imm_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Immunity','FontSize',15)
title('Immunity distribution','FontSize',17)
[p_imm,~,stats_imm] = ranksum(imm_car,imm_non_car,'tail','right');

% Time to first detectable infection: (Fig. 5h, the survival curve in the 
% final version of the manuscript was produced in R, see Simulated_data_analysis.R)
day = 206; % July 25th, ~ day of first follow-up
PCR_pos = zeros(1,N);
event = zeros(1,N); % 0=right censored , 1=event , 2=left censored
thresh = par{1,'lod'}; % use lod of PCR or RDT as the threshold
tmp = zeros(1,N); % 0=PCR- at first follow-up, 1=PCR+ at first follow-up
for i=1:N
    all_par = sum(A{i,4},1);
    ind = find(A{i,3}>=day-par{1,'eds'},1);
    if all_par(ind)>=thresh
        PCR_pos(i) = 0;
        event(i) = 2;
        tmp(i) = 1;
    elseif any(all_par(ind:length(all_par))>thresh) %&& all_par(ind)<thresh
        if treated(i)==1
            if any(all_par(2:size(all_par,2))>thresh)
                PCR_pos(i) = A{i,3}(ind-1+find(all_par(ind:size(all_par,2))>thresh,1))-A{i,3}(ind);
                event(i) = 1;
            else
                PCR_pos(i) = max(A{i,3})-A{i,3}(ind);
                event(i) = 0;
            end
        else
            PCR_pos(i) = A{i,3}(ind-1+find(all_par(ind:length(all_par))>thresh,1))-A{i,3}(ind);
            event(i) = 1;
        end
    else
        PCR_pos(i) = max(A{i,3})-A{i,3}(ind);
        event(i) = 0;
    end
end
% Save survival curve data to visualize survival curves using R (see Simulated_data_analysis.R):
% save('survival-curve-data-1.mat','PCR_pos','day','event','thresh','treated','carrier')
PCR_pos_rp = PCR_pos(logical(treated));
PCR_pos_rnpn = PCR_pos(logical(1-carrier));
PCR_pos_rp = sort(PCR_pos_rp);
PCR_pos_rnpn = sort(PCR_pos_rnpn);
[t_rp,n_rp,m_rp,~,KM_rp] = KM(PCR_pos_rp,365);
[t_rnpn,n_rnpn,m_rnpn,~,KM_rnpn] = KM(PCR_pos_rnpn,365);
stairs(t_rp,KM_rp,'c','LineWidth',2)
hold on
stairs(t_rnpn,KM_rnpn,'r','LineWidth',2)
xlim([0 365])
xlabel('Time [days]','FontSize',15)
ylabel('Fraction without infection','FontSize',15)
title('Time to first detectable infection','FontSize',17)
legend({'Carriers (treated)','Non-carriers'},'Location','northeast','Fontsize',15)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots and analysis of simulation of 1,000 individuals with heterogeneous age and fixed FOI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use "Data-3-x.mat" which contains columns with: age, PMR, t, P, S, C, bite times, and FOI
% (reproduces information in Table S2)

load('Data-3-6.mat')

N = size(A,1); % number of individuals
carrier = zeros(N,1);
treated = zeros(N,1);
for i = 1:N
    carrier(i,1) = sum(A{i,4}(:,1))>par{1,'lod'};
    treated(i,1) = sum(A{i,4}(:,1))>par{1,'lod_rdt'}; % treated if parasites above lod for RDT
end

%%% FOI:
foi = par{1,'bites_mal_max'};

%%% Number of carriers:
n_car = sum(treated);

%%% Time to first bite comparison:
bt_first_all = [];
bt_first_car = [];
bt_first_non_car = [];
for i=1:N
    if ~isempty(A{i,7})
        bt_first_all = [bt_first_all;A{i,7}(1)];
        if treated(i)==1
            bt_first_car = [bt_first_car;A{i,7}(1)];
        end
        if carrier(i)==0
            bt_first_non_car = [bt_first_non_car;A{i,7}(1)];
        end
    end
end
[p_fbt,~,~] = ranksum(bt_first_car,bt_first_non_car);
[p_fbt_car,~,~] = ranksum(bt_first_car,bt_first_non_car,'tail','left');
[p_fbt_non_car,~,~] = ranksum(bt_first_car,bt_first_non_car,'tail','right');

%%% Time to PCR positive comparison:
% TODO

%%% Age comparison:
ages = horzcat(A{:,1})';
age_car = ages(treated==1);
age_non_car = ages(carrier==0);
[p_age,~,~] = ranksum(age_car,age_non_car);
[p_age_non_car,~,~] = ranksum(age_car,age_non_car,'tail','left');
[p_age_car,~,~] = ranksum(age_car,age_non_car,'tail','right');

%%% PMR comparison:
PMRs_all = vertcat(A{:,2})';
PMRs_all = PMRs_all(:);
PMRs_tr = vertcat(A{logical(treated),2});
PMRs_non_car = vertcat(A{logical(1-carrier),2});
[p_pmr,~,~] = ranksum(PMRs_tr,PMRs_non_car);
[p_pmr_car,~,~] = ranksum(PMRs_tr,PMRs_non_car,'tail','left');
[p_pmr_non_car,~,~] = ranksum(PMRs_tr,PMRs_non_car,'tail','right');

%%% Cross-reactive immunity comparison:
imm_all = zeros(N,1);
for i=1:N
    imm_all(i) = A{i,6}(length(A{i,6}));
end
imm_car = imm_all(logical(treated));
imm_non_car = imm_all(logical(1-carrier));
[p_imm,~,~] = ranksum(imm_car,imm_non_car);
[p_imm_car,~,~] = ranksum(imm_car,imm_non_car,'tail','left');
[p_imm_non_car,~,~] = ranksum(imm_car,imm_non_car,'tail','right');

%%% Survival curves: time from first follow-up to PCR+ (Fig. S7)
% (Fig. 5g for Data-3-6.mat)
% The final survival curves were produced using R (see
% Simulated_data_analysis.R).
day = 206; % July 25th, ~ day of first follow-up
PCR_pos = zeros(1,N);
event = zeros(1,N); % 0=right censored , 1=event , 2=left censored
thresh = par{1,'lod'}; % use lod of PCR or RDT as the threshold
tmp = zeros(1,N); % 0=PCR- at first follow-up, 1=PCR+ at first follow-up
for i=1:N
    all_par = sum(A{i,4},1);
    ind = find(A{i,3}>=day-par{1,'eds'},1);
    if all_par(ind)>=thresh
        PCR_pos(i) = 0;
        event(i) = 2;
        tmp(i) = 1;
    elseif any(all_par(ind:length(all_par))>thresh) %&& all_par(ind)<thresh
        if treated(i)==1
            if any(all_par(2:size(all_par,2))>thresh)
                PCR_pos(i) = A{i,3}(ind-1+find(all_par(ind:size(all_par,2))>thresh,1))-A{i,3}(ind);
                event(i) = 1;
            else
                PCR_pos(i) = max(A{i,3})-A{i,3}(ind);
                event(i) = 0;
            end
        else
            PCR_pos(i) = A{i,3}(ind-1+find(all_par(ind:length(all_par))>thresh,1))-A{i,3}(ind);
            event(i) = 1;
        end
    else
        PCR_pos(i) = max(A{i,3})-A{i,3}(ind);
        event(i) = 0;
    end
end
% Save survival curve data to visualize survival curves using R (see Simulated_data_analysis.R):
% save('survival-curve-data-3-1.mat','PCR_pos','day','event','thresh','treated','carrier')
PCR_pos_rp = PCR_pos(logical(treated));
PCR_pos_rnpn = PCR_pos(logical(1-carrier));
PCR_pos_rp = sort(PCR_pos_rp);
PCR_pos_rnpn = sort(PCR_pos_rnpn);
[t_rp,n_rp,m_rp,~,KM_rp] = KM(PCR_pos_rp,365);
[t_rnpn,n_rnpn,m_rnpn,~,KM_rnpn] = KM(PCR_pos_rnpn,365);
stairs(t_rp,KM_rp,'c','LineWidth',2)
hold on
stairs(t_rnpn,KM_rnpn,'r','LineWidth',2)
xlim([0 365])
xlabel('Time [days]','FontSize',15)
ylabel('Probability of being PCR-','FontSize',15)
title('Time from first follow-up to PCR+','FontSize',17)
legend({'Carriers (treated)','Non-carriers'},'Location','northeast','Fontsize',15)
hold off

%%% Boxplot of age distribution (Fig. 5c for Data-3-6.mat)
boxp = boxplot([ages/365;age_car/365;age_non_car/365],...
    [zeros(length(ages),1);ones(length(age_car),1);2*ones(length(age_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Age [years]','FontSize',15)
title('Age distribution','FontSize',17)

%%% Boxplot for cross-reactive immunity (Fig. 5e for Data-3-6.mat)
boxp = boxplot([imm_all;imm_car;imm_non_car],...
    [zeros(length(imm_all),1);ones(length(imm_car),1);2*ones(length(imm_non_car),1)],'Notch','off',...
    'Labels',{'All','Carriers','Non-carriers'},'Colors','k','Symbol','.','OutlierSize',13);
set(gca, 'ActivePositionProperty', 'position','FontSize',15)
set(boxp,'LineWidth', 2);
set(findobj(gcf,'LineStyle','--'),'LineStyle','-') % solid lines
ylabel('Immunity','FontSize',15)
title('Immunity distribution','FontSize',17)


% Fig. 5: homogeneous vs heterogeneous infection risk (age distribution, immunity distribution, time to first 
% detectable infection): 1,000 individuals with FOI 0.024 for homogeneous infection risk



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 3. N individuals with specific biting rates %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate 100,000 individuals with different birthdays in the year.
% 10 different biting rates, 10,000 individuals each.
% No treatment.
% Individuals with positive PCR are classified as "carriers", otherwise
% they are non-carriers.
% Simulate all individuals for 20 years.

for k=1:100 % simulate 100 populations of 1,000 individuals
    for l=1:125 % 125 times 8 individuals
        X = 125;
        r = 14; % initial PMR (per cycle)
        a = 7e-6; % strain specific immunity acquisition rate
        g = 1e-6; % general immunity acquisition rate
        d = 3.7e-4; % general immunity loss rate
        e = 0.056; % initially infected RBCs (in millions)
        lod = 1; % detection threshold (parasites/microlitre)
        lod_rdt = 100; % limit of detection for RDT (parasites/microliter)
        Z_p = 0.005; % 2e-7 is 1 par./5L, parasite threshold (parasites/microlitre)
        liver = 7; % duration of liver-stage (in days)
        eds = 181; % end of dry season (day of the year)
        % March 15th: 74th day      Feburary 15th: 46
        % June 30th:  181st day
        bites_dry = 0; % biting rate in the dry season (in bites per day)
        % simulate 100 individuals for each FOI
        bites_mal_all = [repmat(0.004,X*8/10,1);repmat(0.008,X*8/10,1);repmat(0.012,X*8/10,1);...
            repmat(0.016,X*8/10,1);repmat(0.020,X*8/10,1);repmat(0.024,X*8/10,1);repmat(0.028,X*8/10,1);...
            repmat(0.032,X*8/10,1);repmat(0.036,X*8/10,1);repmat(0.04,X*8/10,1)]; % vector with all biting rates
        
        % table with parameters:
        par = table('Size',[1,13],'VariableTypes',{'double','double','double','double'...
            'double','double','double','double','double','double','double','double','double'},...
            'VariableNames',{'g','r','d','a','e','lod','lod_rdt','Z_p','liver','eds',...
            'bites_dry','bites_mal_min','bites_mal_max'});
        par{1,'g'} = g; par{1,'r'} = r; par{1,'d'} = d; par{1,'a'} = a;
        par{1,'e'} = e; par{1,'lod'} = lod; par{1,'lod_rdt'} = lod_rdt; par{1,'Z_p'} = Z_p;
        par{1,'liver'} = liver; par{1,'eds'} = eds; par{1,'bites_dry'} = bites_dry;
        
        % blood volume (in litres) depending on age (in days)
        V =@(a) (a<22*365).*(0.00059.*a+0.3) + (a>=22*365).*5;
        
        
        N = 8; % number of individuals
        n = 500;
        age_max = 20*365; % surveillance period (days)
        
        A_tmp = cell(N,1); % save birth, bites_mal_ind, t, P, S, G, bt for each individual
        
        parfor i=1:N
            birth = rand(1)*365+1; % time of birth (for bite rate)
            if birth==366
                birth = 1;
            end
            bites_mal_ind = bites_mal_all((l-1)*N+i); % biting rate during the malaria
            
            bt = zeros(1,100*age_max/365); % bite times
            t = zeros(3100*age_max/365,1); % time in days
            P = zeros(n,3100*age_max/365); % parasite concentration
            S = zeros(n,3100*age_max/365); % strain specific immunity
            C = zeros(1,3100*age_max/365); % general immunity
            
            % first bite:
            tcb = ttnb(birth,bites_dry,bites_mal_ind); % time to first bite
            tnb = tcb + ttnb(birth+tcb,bites_dry,bites_mal_ind); % time of next bite
            bt(1) = tcb; % bite times
            
            y0 = zeros(2*n+1,1);
            y0(1,1) = e/V(tcb);  % transfer of parasites by bite
            [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+liver) (tnb+liver)],y0);
            t(2:length(tmp),1) = tmp(1:(length(tmp)-1));
            P(:,2:length(tmp)) = (y(1:(length(tmp)-1),1:n))';
            S(:,2:length(tmp)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
            C(1,2:length(tmp)) = (y(1:(length(tmp)-1),2*n+1))';
            tcb = tnb;
            ind = 1;
            ind2 = length(tmp)+1;
            
            % later bites
            while(tcb+liver<age_max)
                ind = ind+1;
                bt(ind) = tcb; % save time of current bite
                tnb = tcb + ttnb(birth+tcb,bites_dry,bites_mal_ind); % time of next bite
                
                % if parasites are under a threshold, set parasite concentration to 0
                % also set strain specific immunity to zero
                % (it would only decay exponentially and not influence G or other parasite strains)
                if any(P(1:(find(max(P,[],2)==0,1)-1),find(t==bt(ind-1)+liver):(ind2-1))<Z_p,'all')
                    for j=1:(find(max(P,[],2)==0,1)-1)
                        tmp = find(P(j,find(t==bt(ind-1)+liver):(ind2-1))<Z_p,1);
                        P(j,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
                        S(j,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
                    end
                end
                
                % initializing:
                y0 = [P(:,ind2-1);S(:,ind2-1);C(1,ind2-1)];
                
                % inoculate parasites in first "empty" row
                cs =1;
                if any(P(:,ind2-1)==0)
                    cs = find(P(:,ind2-1)==0,1);
                else
                    break
                end
                y0(cs,1) = y0(cs,1) + e/V(tcb);
                
                if tnb+liver <= age_max
                    [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+liver) (tnb+liver)],y0);
                    t(ind2:(ind2+length(tmp)-2),1) = tmp(1:(length(tmp)-1));
                    P(:,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),1:n))';
                    S(:,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),(n+1):(2*n)))';
                    C(1,ind2:(ind2+length(tmp)-2)) = (y(1:(length(tmp)-1),2*n+1))';
                    tcb = tnb;
                    ind2 = ind2+length(tmp)-1;
                else
                    [tmp,y] = ode45(@(t,y) intrahost(y(1:n,1),y((n+1):(2*n),1),y(2*n+1,1),par),[(tcb+liver) age_max],y0);
                    t(ind2:(ind2+length(tmp)-1),1) = tmp(1:(length(tmp)));
                    P(:,ind2:(ind2+length(tmp)-1)) = (y(1:(length(tmp)),1:n))';
                    S(:,ind2:(ind2+length(tmp)-1)) = (y(1:(length(tmp)),(n+1):(2*n)))';
                    C(1,ind2:(ind2+length(tmp)-1)) = (y(1:(length(tmp)),2*n+1))';
                    tcb = tnb;
                    ind2 = ind2+length(tmp);
                end
            end
            
            % if parasites are under a threshold, set parasite concentration to 0
            % also set strain specific immunity to zero
            % (it would only decay exponentially and not influence G or other parasite strains)
            if any(P(1:(find(max(P,[],2)==0,1)-1),find(t==bt(ind-1)+liver):(ind2-1))<Z_p,'all')
                for j=1:(find(max(P,[],2)==0,1)-1)
                    tmp = find(P(j,find(t==bt(ind-1)+liver):(ind2-1))<Z_p,1);
                    P(j,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
                    S(j,(find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)) = zeros(size((find(t==bt(ind-1)+liver)+tmp-1):(ind2-1)));
                end
            end
            
            rows = find(max(P,[],2)==0,1)-1;
            bt = bt(1:ind);
            t = t(1:ind2-1,1);
            P = P(1:rows,1:(ind2-1));
            S = S(1:rows,1:(ind2-1));
            C = C(1,1:(ind2-1));
            
            A_tmp{i,1}={birth,t,P,S,C,bt,bites_mal_ind};
            
            % progress:
            %         disp(i)
        end
        
        A = cell(N,7);
        for i=1:N
            for j=1:7
                A{i,j} = A_tmp{i}{j};
            end
        end
        
        filename = ['Data-2-prel-' num2str(l) '.mat'];
        save(filename,'A','par')
        
        % progress:
        %     disp('.')
        
        clearvars('-except', 'k');
        
    end
    
    % put together the l files into one file:
    A_tmp = cell(104,7);
    
    X=125;
    for l=1:X
        filename = ['Data-2-prel-' num2str(l) '.mat'];
        load(filename)
        
        for i=1:8
            for j=1:7
                A_tmp{i+(l-1)*8,j} = A{i,j};
            end
        end
        
    end
    
    A = A_tmp;
    filename = ['Data-2-' num2str(k) '.mat'];
    save(filename,'A','par')
    
    clearvars('-except', 'k');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots and analysis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
N = 100;
A_all = cell(100000,7);
for k=1:N
    filename = ['Data-2-' num2str(k) '.mat'];
    load(filename)
    for i=1:1000
        % columns in A: birth time, time, parasites, strain specific immunity,
        % general immunity, bite times, FOI
        for j=[1:3,7] % only these columns are needed for the plots
            A_all{((k-1)*1000+i),j} = A{i,j};
        end
    end
    disp(".")
    if mod(k,10)==0
        disp(k)
    end
end
A = A_all;

% load parameters:
g = par{1,'g'}; r = par{1,'r'}; d = par{1,'d'}; a = par{1,'a'};
e = par{1,'e'}; lod = par{1,'lod'}; lod_rdt = par{1,'lod_rdt'}; Z_p = par{1,'Z_p'};
liver = par{1,'liver'}; eds = par{1,'eds'}; bites_dry = par{1,'bites_dry'};

N = size(A,1);
bites_mal_all = vertcat(A{:,7});
age_max = 20;

%%% Age of first carriage boxplot (Fig. S4)
age_first_car = repmat(age_max+1,N,1);
for i=1:N
    birth_tmp = A{i,1};
    t_tmp = A{i,2};
    P_tmp = A{i,3};
    
    x = mod(eds-birth_tmp-1,365)+1; % age at the last day of the dry season
    for j=0:floor((max(t_tmp)-x)/365)
        [~,ind] = min(abs(t_tmp-(x+j*365)));
        if sum(P_tmp(:,ind))>lod
            age_first_car(i) = t_tmp(ind)/365;
            %            age_first_car(i) = floor(t_tmp(ind)/365);
            break
        end
    end
end

boxplot(age_first_car,bites_mal_all)
xlabel('Force Of Infection (FOI) [bites/day]')
ylabel('Age of first parasite carriage [years]')
title('Age of first parasite carriage at the end of the dry season')


%%% Age of first carriage heatmap (Fig. S3)
N = size(A,1);
bites_mal_all = vertcat(A{:,7});
age_max = 20;
first_carriage = zeros(length(unique(bites_mal_all)),age_max+1);
ages = 0:age_max;
foi = sort(unique(bites_mal_all),'descend');

age_first_car = repmat(age_max,N,1);
for i=1:N
    birth_tmp = A{i,1};
    t_tmp = A{i,2};
    P_tmp = A{i,3};
    
    x = mod(eds-birth_tmp-1,365)+1; % age at the last day of the dry season
    for j=0:floor((max(t_tmp)-x)/365)
        [~,ind] = min(abs(t_tmp-(x+j*365)));
        if sum(P_tmp(:,ind))>lod
            %            age_first_car(i) = t_tmp(ind)/365;
            age_first_car(i) = floor(t_tmp(ind)/365);
            break
        end
    end
end

for i=1:length(unique(bites_mal_all))
    for j=1:(age_max+1)
        first_carriage(i,j)=sum(age_first_car==ages(j) & bites_mal_all==foi(i))/(N/10);
    end
end
first_carriage = round(first_carriage*10000)/100;

heatmap(ages,foi,first_carriage)
xlabel('Age of first parasite carriage [years]')
ylabel('Force Of Infection (FOI) [bites/day]')
title('Age of first parasite carriage at the end of the dry season')
% colormap(parula)


%%% Fraction of carriers by age and FOI: (Fig. 3)
N = size(A,1);
eds = 181; % end of dry season (day of the year)
lod = 1; % PCR detection threshold (parasites/microlitre)
age_groups = {'0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10',...
    '10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18',...
    '18-19','19-20'};
age_all = zeros(N*length(age_groups),1);
age_group_all = strings(N*length(age_groups),1);
foi_all = zeros(N*length(age_groups),1);
PCR_pos_all = zeros(N*length(age_groups),1);

for i=1:N % for each individual
    birth_tmp = A{i,1};
    t_tmp = A{i,2};
    P_tmp = A{i,3};
    
    % age at the last day of the first experienced dry season
    x = mod(eds-birth_tmp,365);
    
    heatmap_tmp = A{i,7}; % FOI for individual i
    for j = 1:length(age_groups) % for each age
        age_tmp = x+(j-1)*365; % age at the end of the jth dry season
        age_group_tmp = age_groups{floor(age_tmp/365)+1}; % age group
        [~,ind] = min(abs(t_tmp-age_tmp));
        PCR_pos_tmp = sum(P_tmp(:,ind),1)>lod;
        
        age_all((i-1)*length(age_groups)+j) = age_tmp;
        age_group_all((i-1)*length(age_groups)+j) = age_group_tmp;
        foi_all((i-1)*length(age_groups)+j) = heatmap_tmp;
        PCR_pos_all((i-1)*length(age_groups)+j) = PCR_pos_tmp;
    end
    
    if mod(i,1000)==0
        disp(i)
    end
end

tbl = table(age_group_all,foi_all,PCR_pos_all);

foi = sort(unique(foi_all),'descend');
age_groups = {'0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10',...
    '10-11','11-12','12-13','13-14','14-15','15-16','16-17','17-18',...
    '18-19','19-20'};
heatmap_tmp = zeros(length(unique(foi_all)),length(unique(age_group_all)));
for i=1:length(unique(foi_all)) % for each foi
    for j=1:length(unique(age_group_all)) % for each age group
        heatmap_tmp(i,j) = sum(PCR_pos_all(foi_all==foi(i) & age_group_all==age_groups(j)))/10000;
    end
end

mycolors = [sscanf('E6F1F8','%2x%2x%2x',[1 3])/255; sscanf('CCE3F2','%2x%2x%2x',[1 3])/255;...
    sscanf('B3D5EB','%2x%2x%2x',[1 3])/255; sscanf('9DC9E6','%2x%2x%2x',[1 3])/255; ...
    sscanf('83BBDF','%2x%2x%2x',[1 3])/255;sscanf('6AACD8','%2x%2x%2x',[1 3])/255;...
    sscanf('509ED2','%2x%2x%2x',[1 3])/255;sscanf('3790CB','%2x%2x%2x',[1 3])/255;...
    sscanf('1D82C5','%2x%2x%2x',[1 3])/255;sscanf('0474BE','%2x%2x%2x',[1 3])/255];
colormap (mycolors);
contourf(0.5:1:19.5,flip(0.004:0.004:0.04),heatmap_tmp,[0:0.1:1])%,'ShowText','on')
colorbar
xlabel('Age [years]')
ylabel('Force Of Infection (FOI) [bites/day]')
yticks(sort(unique(foi_all)))
title('Fraction of carriers')
xlim([0.5,19.5])
ylim([min(foi_all),max(foi_all)])
hold off


%%% Number of bites in the previous transmission season for carriers and
% non-carriers:
% load data
N = 100;
A_all = cell(100000,7);
for k=1:N
    filename = ['Data-2-' num2str(k) '.mat'];
    load(filename)
    for i=1:1000
        % columns in A: birth time, time, parasites, strain specific immunity,
        % general immunity, bite times, FOI
        for j=[1:3,6,7] % only these columns are needed for the plots
            A_all{((k-1)*1000+i),j} = A{i,j};
        end
    end
    disp(".")
    if mod(k,10)==0
        disp(k)
    end
end
A = A_all;

% load parameters:
g = par{1,'g'}; r = par{1,'r'}; d = par{1,'d'}; a = par{1,'a'};
e = par{1,'e'}; lod = par{1,'lod'}; lod_rdt = par{1,'lod_rdt'}; Z_p = par{1,'Z_p'};
liver = par{1,'liver'}; eds = par{1,'eds'}; bites_dry = par{1,'bites_dry'};

N = size(A,1);
eds = 181; % end of dry season (day of the year)
bites_mal_all = vertcat(A{:,7});
age_max = 20;

% age structure: uniform distribution between 0 and 20 years
ages = unidrnd(age_max,N,1)-1;
carrier = zeros(N,1); % 1 for carriers, 0 otherwise
num_bites = zeros(N,1); % number of infectious bites in the previous season
last_bite_time = nan*zeros(N,1); % time between last bite and the end of the transmission season

for i=1:100000
    age_tmp = ages(i);
    birth_tmp = A{i,1};
    t_tmp = A{i,2};
    P_tmp = A{i,3};
    bt_tmp = A{i,6};
    % find the time of the last day of the dry season with the given age:
    x = mod(eds-birth_tmp-1,365)+1+365*age_tmp;
    if t_tmp(2)<x % experienced at least one bite
        [~,ind_tmp] = min(abs(t_tmp(t_tmp<x)-x)); % last time point in the dry season
    else
        ind_tmp = 1;
    end
    % carriage:
    if sum(P_tmp(:,ind_tmp))>lod
        carrier(i,1) = 1;
    end
    % number of infectious bites in the previous season:
    num_bites(i,1) = sum(bt_tmp<t_tmp(ind_tmp) & bt_tmp>t_tmp(ind_tmp)-365);
    if x<eds % never experienced a transmission season
        num_bites(i,1) = nan;
    end
    % time between last bite and the end of the transmission season
    x_tr = x-181; % last day of the previous transmission season
    if any(bt_tmp<x_tr)
        last_bite_time(i,1) = x_tr-max(bt_tmp(bt_tmp<x_tr));
    end
    
    if mod(i,1000)==0
        disp(i)
    end
end

data = [bites_mal_all,ages,carrier,num_bites,last_bite_time];
save('Data-2-car-vs-non-car','data')

% load data
load('Data-2-car-vs-non-car.mat')

% This data is analyzed in R using the script "Simulated_data_analysis.R"


