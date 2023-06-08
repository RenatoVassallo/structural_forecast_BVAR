%% MATLAB CODE FOR:
% ========================================================================  
%                FORECASTING AND NOWCASTING WITH TEXT AS DATA
%                     Economic Policy Uncertainty Shocks 
%                           
%                              June, 2023
%                           
% ========================================================================
% AUTHOR:  Renato Vassallo
% CONTACT: renato.vassallo@bse.eu
%-------------------------------------------------------------------------
clear; clc;
addpath('functions');
addpath('data');

% ==============================
%% STEP 0: DATA PRE-PROCESSING
% ==============================
% Load US data
data = readtable('data\homework.xlsx');
data.Dates = datetime(data.Dates);
data.Three_Component_Index = [];

% Rename columns
data.Properties.VariableNames{'PAYEMS'}   = 'EMP';
data.Properties.VariableNames{'INDPRO'}   = 'IND';
data.Properties.VariableNames{'FEDFUNDS'} = 'FED';
data.Properties.VariableNames{'News_Based_Policy_Uncert_Index'} = 'EPU';

% Compute MoM growth rates
data.PCE_Growth  = 100*[0; diff(data.PCE)] ./ data.PCE;
data.EMP_Growth  = 100*[0; diff(data.EMP)] ./ data.EMP;
data.IND_Growth  = 100*[0; diff(data.IND)] ./ data.IND;

% Filter the data for the specified period
start_date = datetime('01-Feb-1985');
end_date = datetime('01-Dec-2019');
filtered_data = data(data.Dates >= start_date & data.Dates <= end_date, :);

% Extract the variables
EPU = filtered_data.EPU;
IND_Growth = filtered_data.IND_Growth;
PCE_Growth = filtered_data.PCE_Growth;
EMP_Growth = filtered_data.EMP_Growth;

% Create the final dataset
d = [EPU, IND_Growth, PCE_Growth, EMP_Growth];

N = cols(d);
L = 1; %lag length of the VAR
Y = d;
%take lags
X = [];

for j = 1:L
    X = [X lag0(d,j) ];
end

X = [ones(rows(X),1) X];
Y = Y(L+1:end,:);
X = X(L+1:end,:);
T = rows(X);

nhor = 20; %IRF horizon
REPS = 40000;
BURN = 30000;

% ==========================
%% STEP 1: MINNESOTA PRIOR
% ==========================

% Set priors and starting values
%compute standard deviation of each series residual via an ols regression
%to be used in setting the prior
std = zeros(N, 1);
for i = 1:N
    y = Y(:, i);
    x = X(:, [1, i+1]);
    
    b0 = inv(x' * x) * (x' * y);
    residuals = y - x * b0;
    s = sqrt((residuals' * residuals) / (length(y) - 2));
    std(i) = s;
end

%parameters to control the prior
lamda1 = 0.2;  %tightness prior on own lags
lamda2 = 0.5;  %tightness prior on lags of variables other than the dependent variable
lamda3 = 1;    %tightness of prior on higher lags 
lamda4 = 10^5; %tightness of prior on the constant term

%specify the prior mean of the coefficients of the Two equations of the VAR
B0=zeros((N*L+1),N);
for i=1:N
    B0(i+1,i) = 0.75;
end
B0=vec(B0);

%Specify the prior variance of vec(B)
H=zeros(20,20);
%for EPU equation of the VAR
H(1,1)=(std(1)*lamda4)^2;                  %constant
H(2,2)=(lamda1)^2;                         %own lag
H(3,3)=((std(1)*lamda1*lamda2)/std(2))^2;  %lag of ind
H(4,4)=((std(1)*lamda1*lamda2)/std(3))^2;  %lag of inf
H(5,5)=((std(1)*lamda1*lamda2)/std(4))^2;  %lag of emp
%for INDPRO equation of the VAR
H(6,6)=(std(2)*lamda4)^2;                %constant
H(7,7)=((std(2)*lamda1*lamda2)/std(1))^2;  %lag of epu
H(8,8)=(lamda1)^2;                         %own lag
H(9,9)=((std(2)*lamda1*lamda2)/std(3))^2;  %lag of inf
H(10,10)=((std(2)*lamda1*lamda2)/std(4))^2;  %lag of emp
%for INFL equation of the VAR
H(11,11)=(std(3)*lamda4)^2;                  %constant
H(12,12)=((std(3)*lamda1*lamda2)/std(1))^2;  %lag of epu
H(13,13)=((std(3)*lamda1*lamda2)/std(2))^2;  %lag of ind
H(14,14)=(lamda1)^2;                         %own lag
H(15,15)=((std(3)*lamda1*lamda2)/std(4))^2;  %lag of emp
%for EMPL equation of the VAR
H(16,16)=(std(4)*lamda4)^2;                  %constant
H(17,17)=((std(4)*lamda1*lamda2)/std(1))^2;  %lag of epu
H(18,18)=((std(4)*lamda1*lamda2)/std(2))^2;  %lag of ind
H(19,19)=((std(4)*lamda1*lamda2)/std(3))^2;  %lag of emp
H(20,20)=(lamda1)^2;                         %own lag

%prior scale matrix for sigma the VAR covariance
S=eye(N);
%prior degrees of freedom
alpha=N+1;

%starting values for the Gibbs sampling algorithm
Sigma=eye(N);
betaols=vec(inv(X'*X)*(X'*Y));


% ==================================
%% STEP 2: Start the Gibbs Sampler
% ==================================

coef1=[]; %will store coeff for VAR parameters
coef2=[]; %will store variances
out1=[];  %will store IRFs for EPU
out2=[];  %will store IRFs for IP
out3=[];  %will store IRFs for PCE
out4=[];  %will store IRFs for EMP

for j=1:REPS

%draw the VAR coefficients
M=inv(inv(H)+kron(inv(Sigma),X'*X))*(inv(H)*B0+kron(inv(Sigma),X'*X)*betaols);
V=inv(inv(H)+kron(inv(Sigma),X'*X));

%check for stability of the VAR
check=-1;
while check<0
beta=M+(randn(1,N*(N*L+1))*chol(V))';
CH=stability(beta,N,L);
if CH==0
    check=10;
end
end

%draw sigma from the IW distribution
e=Y-X*reshape(beta,N*L+1,N);
%scale matrix
scale=e'*e+S;
Sigma=IWPQ(T+alpha,inv(scale));

if j>BURN
    %impulse response using a cholesky decomposition
    A0=chol(Sigma);
    v=zeros(nhor,N);
    v(L+1,1) = 3; %shock the EPU

    yhat=zeros(nhor,N);
    for i=L+1:nhor
        yhat(i,:)=[0 yhat(i-1,:)]*reshape(beta,N*L+1,N)+v(i,:)*A0;
    end

    coef1=[coef1;beta'];
    coef2=[coef2;Sigma];
    out1=[out1 yhat(L+1:end,1)];
    out2=[out2 yhat(L+1:end,2)];
    out3=[out3 yhat(L+1:end,3)];
    out4=[out4 yhat(L+1:end,4)];

end
end

% ==========================
%% STEP 3: GENERATE GRAPHS
% ==========================


% FIGURE 1: VARIABLES OF THE MODEL
%----------------------------------
x0=100;
y0=100;
width=1000; %largo
height=650; %alto

figure('Name','Variables of the Model','NumberTitle','off')
t = tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact');

titles = {'EPU (Index)', 'Industrial Production (M-o-M growth)', 'PCE (M-o-M growth)', 'Employment (M-o-M growth)'};
variables = {'EPU', 'IND_Growth', 'PCE_Growth', 'EMP_Growth'};

for i = 1:N
    nexttile
    plot(filtered_data.Dates, filtered_data.(variables{i}));
    hold on;
    line([filtered_data.Dates(1), filtered_data.Dates(end)], [0, 0], 'Color', 'k', 'LineStyle', ':');
    hold off; axis tight;
    title(titles{i});
end
set(gcf,'position',[x0,y0,width,height])
exportgraphics(t, 'variables.png', 'BackgroundColor','white','Resolution', 300);


% FIGURE 2: MARGINAL POSTERIOR DISTRIBUTIONS
%--------------------------------------------
x0=100;
y0=100;
width=1000; %largo
height=650; %alto

figure('Name','Marginal Posterior Distributions','NumberTitle','off')
t = tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact');

titles = {'$\phi_{11}$', '$\phi_{21}$', '$\phi_{31}$', '$\phi_{41}$'};
positions = [2, 7, 12, 17];
histColor = [0.5 0.5 0.5]; % Specify the desired color for the histogram

for i = 1:4
    nexttile
    data = coef1(:, positions(i));
    histogram(data,'FaceColor',[0.0 0.45 0.74]);
    axis tight;
    title(titles(i),'Interpreter', 'latex','FontSize',14);

    medianLine = median(data);
    hold on;
    line([medianLine medianLine], ylim, 'Color', 'r', 'LineWidth', 1);
    hold off;
end
set(gcf,'position',[x0,y0,width,height])
exportgraphics(t, 'coef_distributions.png', 'BackgroundColor','white','Resolution', 300);



% FIGURE 3: IMPULSE-RESPONSE FUNCTIONS TO AN EPU SHOCK
%------------------------------------------------------

x0=100;
y0=100;
width=1000; %largo
height=650; %alto

figure('Name','Impulse Responses to a 1 s.d shock to EPU equation','NumberTitle','off')
t = tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact');

data = {out1, out2, out3, out4};
titles = {'Response of EPU', 'Response of Industrial Production Growth', 'Response of Consumption Growth', 'Response of Employment Growth'};
percentiles = [50, 16, 84]; % Percentiles for filling

for i = 1:4
    nexttile
    prc = prctile(data{i}, percentiles, 2);
    x = [1:size(prc, 1), fliplr(1:size(prc, 1))]; % x-values for filling
    y = [prc(:, 2)', fliplr(prc(:, 3)')]; % y-values for filling
    grid on;
    hold on;
    fill(x, y, [1, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Fill area between percentiles
    plot(prc(:, 1), 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5);
    plot(zeros(size(prc, 1), 1), 'k','LineStyle','-','LineWidth', 0.25);
    hold off;
    title(titles{i},"FontSize",12);
    axis tight;
end
leg = legend({'68% Credible Regions', 'Median Response'},'Orientation', 'horizontal','FontSize', 12);
leg.Layout.Tile = 'south';
set(gcf,'position',[x0,y0,width,height])
exportgraphics(t, 'irf_epu.png', 'BackgroundColor','white','Resolution', 300);


    