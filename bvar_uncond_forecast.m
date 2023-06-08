%% MATLAB CODE FOR:
% ========================================================================  
%                FORECASTING AND NOWCASTING WITH TEXT AS DATA
%                          Unconditional Forecast 
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

% Data till 2023-03-01
EPU_all = data.EPU;
IND_Growth_all = data.IND_Growth;
PCE_Growth_all = data.PCE_Growth;
EMP_Growth_all = data.EMP_Growth;

% Create the final dataset
d = [EPU, IND_Growth, PCE_Growth, EMP_Growth];
d_all = [EPU_all, IND_Growth_all, PCE_Growth_all, EMP_Growth_all];
N = cols(d);
N_all = cols(d_all);

L = 1; %lag length of the VAR
Y = d;
Y_all = d_all;

%take lags
X=[];
for j=1:L
    X=[X lag0(d,j)];
end

X = [ones(rows(X),1) X];
Y = Y(L+1:end,:);
Y_all = Y_all(L+1:end,:);
X = X(L+1:end,:);
T = rows(X);

horizon = 12; %forecast horizon

%constrained values for EPU
path = [306;328;327;261;188;150;150;150;150;150;150;150]; %ones(12,1).*160; 

% =================================
%% STEP 1: MODEL ESTIMATION
% =================================

B = X\Y;  %ols estimate
res = Y-X*B;
sigma = (res'*res)/T;
A0 = chol(sigma); 

%calculate impulse responses to be used to construct R
S=zeros(1,N);
S(1)=1; %shock to 1st eq
Z1=irfsim(B,N,L,A0,S,horizon+L);
S=zeros(1,N);
S(2)=1; %shock to 2nd eq
Z2=irfsim(B,N,L,A0,S,horizon+L);
S=zeros(1,N);
S(3)=1; %shock to 3rd eq
Z3=irfsim(B,N,L,A0,S,horizon+L);
S=zeros(1,N);
S(4)=1; %shock to 4th eq
Z4=irfsim(B,N,L,A0,S,horizon+L);

%calculate unconditional forecast to be used to construct r
yhat1=zeros(horizon+L,N);
yhat1(1:L,:)=Y_all(end-L+1:end,:);
for i=L+1:horizon+L   
    x=[];
    for j=1:L
    x=[x yhat1(i-j,:)];
    end
    yhat1(i,:)=[1 x]*B;  
end
yhat1=yhat1(L+1:end,:);

%construct the R matrix
ZT1  = [Z1(1,1)  Z2(1,1)  Z3(1,1)  Z4(1,1)];
ZT2  = [Z1(2,1)  Z2(2,1)  Z3(2,1)  Z4(2,1)];
ZT3  = [Z1(3,1)  Z2(3,1)  Z3(3,1)  Z4(3,1)];
ZT4  = [Z1(4,1)  Z2(4,1)  Z3(4,1)  Z4(4,1)];
ZT5  = [Z1(5,1)  Z2(5,1)  Z3(5,1)  Z4(5,1)];
ZT6  = [Z1(6,1)  Z2(6,1)  Z3(6,1)  Z4(6,1)];
ZT7  = [Z1(7,1)  Z2(7,1)  Z3(7,1)  Z4(7,1)];
ZT8  = [Z1(8,1)  Z2(8,1)  Z3(8,1)  Z4(8,1)];
ZT9  = [Z1(9,1)  Z2(9,1)  Z3(9,1)  Z4(9,1)];
ZT10 = [Z1(10,1) Z2(10,1) Z3(10,1) Z4(10,1)];
ZT11 = [Z1(11,1) Z2(11,1) Z3(11,1) Z4(11,1)];
ZT12 = [Z1(12,1) Z2(12,1) Z3(12,1) Z4(12,1)];

R=[ZT1 zeros(1,N*horizon-N);
   ZT2 ZT1 zeros(1,N*horizon-2*N);
   ZT3 ZT2 ZT1 zeros(1,N*horizon-3*N);
   ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-4*N);
   ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-5*N);
   ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-6*N);
   ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-7*N);
   ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-8*N);
   ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-9*N);
   ZT10 ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-10*N);
   ZT11 ZT10 ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-11*N);
   ZT12 ZT11 ZT10 ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 ];

%construct the r matrix
r=path-yhat1(:,1);
%compute the restricted structural shocks
ehat=R'*pinv(R*R')*r;
ehat=reshape(ehat,N,horizon)';
%compute the conditional forecast
yhat2=zeros(horizon+L,N);
yhat2(1:L,:)=Y_all(end-L+1:end,:);
for i=L+1:horizon+L   
    x=[];
    for j=1:L
    x=[x yhat2(i-j,:)];
    end
    yhat2(i,:)=[1 x]*B+ehat(i-L,:)*A0;  
end
yhat2=yhat2(L+1:end,:);

% =================================
%% STEP 2: GIBBS SAMPLING ALGORITHM
% =================================
REPS = 15000;
BURN = 5000;

out1unc=[]; %will hold unc forecast for EPU
out2unc=[]; %will hold unc forecast for INDPRO
out3unc=[]; %will hold unc forecast for PCE
out4unc=[]; %will hold unc forecast for EMP
out1con=[]; %will hold cond forecast for EPU
out2con=[]; %will hold cond forecast for INDPRO
out3con=[]; %will hold cond forecast for PCE
out4con=[]; %will hold cond forecast for EMP
yhatg=yhat2; %initialise conditional forecast
sig=sigma; %initialise error covariance
for igibbs=1:REPS
    
    %step 1 DRAW VAR parameters
    datag=[d;yhatg]; %appended data
    YSTAR=datag;
    %take lags
    XSTAR=[];
    for j=1:L
        XSTAR=[XSTAR lag0(datag,j)];
    end
    XSTAR=[ones(rows(XSTAR),1) XSTAR];
    YSTAR=YSTAR(L+1:end,:);
    XSTAR=XSTAR(L+1:end,:);
    T=rows(XSTAR);
    %conditional mean
    M=vec(XSTAR\YSTAR);
    %conditional variance
    V=kron(sig,inv(XSTAR'*XSTAR));
    bg=M+(randn(1,N*(N*L+1))*chol(V))';
    bg1=reshape(bg,N*L+1,N);
    %draw sigma from the IW distribution
    e=YSTAR-XSTAR*bg1;
    scale=e'*e;
    sig=IWPQ(T,inv(scale));
    %A0 matrix
    A0g=chol(sig);


    %step 2 Construct conditional forecast 

    %impulse responses
    S=zeros(1,N);
    S(1)=1; %shock to 1st eq
    Z1=irfsim(bg1,N,L,A0g,S,horizon+L);
    S=zeros(1,N);
    S(2)=1; %shock to 2nd eq
    Z2=irfsim(bg1,N,L,A0g,S,horizon+L);
    S=zeros(1,N);
    S(3)=1; %shock to 3rd eq
    Z3=irfsim(bg1,N,L,A0g,S,horizon+L);
    S=zeros(1,N);
    S(4)=1; %shock to 4th eq
    Z4=irfsim(bg1,N,L,A0g,S,horizon+L);

    %calculate unconditional forecast to be used to construct r
    yhat1=zeros(horizon+L,N);
    yhat1(1:L,:)=Y_all(end-L+1:end,:);
    for i=L+1:horizon+L   
        x=[];
        for j=1:L
        x=[x yhat1(i-j,:)];
        end
        yhat1(i,:)=[1 x]*bg1;  
    end
    yhat1=yhat1(L+1:end,:);
    
    %construct the R matrix
    ZT1  = [Z1(1,1)  Z2(1,1)  Z3(1,1)  Z4(1,1)];
    ZT2  = [Z1(2,1)  Z2(2,1)  Z3(2,1)  Z4(2,1)];
    ZT3  = [Z1(3,1)  Z2(3,1)  Z3(3,1)  Z4(3,1)];
    ZT4  = [Z1(4,1)  Z2(4,1)  Z3(4,1)  Z4(4,1)];
    ZT5  = [Z1(5,1)  Z2(5,1)  Z3(5,1)  Z4(5,1)];
    ZT6  = [Z1(6,1)  Z2(6,1)  Z3(6,1)  Z4(6,1)];
    ZT7  = [Z1(7,1)  Z2(7,1)  Z3(7,1)  Z4(7,1)];
    ZT8  = [Z1(8,1)  Z2(8,1)  Z3(8,1)  Z4(8,1)];
    ZT9  = [Z1(9,1)  Z2(9,1)  Z3(9,1)  Z4(9,1)];
    ZT10 = [Z1(10,1) Z2(10,1) Z3(10,1) Z4(10,1)];
    ZT11 = [Z1(11,1) Z2(11,1) Z3(11,1) Z4(11,1)];
    ZT12 = [Z1(12,1) Z2(12,1) Z3(12,1) Z4(12,1)];

    R=[ZT1 zeros(1,N*horizon-N);
       ZT2 ZT1 zeros(1,N*horizon-2*N);
       ZT3 ZT2 ZT1 zeros(1,N*horizon-3*N);
       ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-4*N);
       ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-5*N);
       ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-6*N);
       ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-7*N);
       ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-8*N);
       ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-9*N);
       ZT10 ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-10*N);
       ZT11 ZT10 ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 zeros(1,N*horizon-11*N);
       ZT12 ZT11 ZT10 ZT9 ZT8 ZT7 ZT6 ZT5 ZT4 ZT3 ZT2 ZT1 ];
    
    %construct the r matrix
    r=path-yhat1(:,1);
    
    %compute the mean of the distribution of restricted structural shocks
    MBAR=R'*pinv(R*R')*r;
    %compute the variance of the distribution of restricted structural shocks
    VBAR=R'*pinv(R*R')*R;
    VBAR=eye(cols(VBAR))-VBAR;
    %draw structural shocks from the N(MBAR,VBAR) distribution
    edraw=MBAR+(randn(1,rows(MBAR))*real(sqrtm(VBAR)))';
    edraw=reshape(edraw,N,horizon)';
    %conditional forecast using new draw of shocks
    yhatg=zeros(horizon+L,N);
    yhatg(1:L,:)=Y_all(end-L+1:end,:);
    for i=L+1:horizon+L   
        x=[];
        for j=1:L
        x=[x yhatg(i-j,:)];
        end
        yhatg(i,:)=[1 x]*bg1+edraw(i-L,:)*A0g;  
    end
    yhatg=yhatg(L+1:end,:);
    
    if igibbs>BURN
        out1unc = [out1unc;[Y_all(:,1);yhat1(:,1)]'];
        out2unc = [out2unc;[Y_all(:,2);yhat1(:,2)]'];
        out3unc = [out3unc;[Y_all(:,3);yhat1(:,3)]'];
        out4unc = [out4unc;[Y_all(:,4);yhat1(:,4)]'];
        out1con = [out1con;[Y_all(:,1);yhatg(:,1)]'];
        out2con = [out2con;[Y_all(:,2);yhatg(:,2)]'];
        out3con = [out3con;[Y_all(:,3);yhatg(:,3)]'];
        out4con = [out4con;[Y_all(:,4);yhatg(:,4)]'];
    end
end

outunc = cat(3,out1unc, out2unc, out3unc, out4unc);
outcon = cat(3,out1con, out2con, out3con, out4con);



% ========================
%% STEP 3: GENERATE GRAPHS
% ========================

% FIGURE 1: UNCONDITIONAL FORECAST (FAN CHART)
%----------------------------------------------

x0=100;
y0=100;
width=1200; %largo
height=700; %alto

figure('Name','Unconditional Forecast','NumberTitle','off')
t = tiledlayout(2,2,'Padding', 'none', 'TileSpacing', 'compact');

titles = {'EPU Index','Industrial Production Growth','Consumption Growth','Employment Growth'};
TT = datetime(2021,06,01):calmonths(1):datetime(2024,03,01);
percentiles = [10 20 30 40 50 60 70 80 90];
fan_color = [0.39 0.2 1];

% Compute the percentiles for each matrix in 'outcon' and 'outunc'
pct_values_outunc = cell(1, size(outunc, 3));
for i = 1:size(outunc, 3)
    pct_values_outunc{i} = prctile(outunc(:, end-size(TT,2)+1:end, i), percentiles);
end

% Loop over each matrix in 'outcon' and 'outunc' and create the corresponding subplot
for i = 1:size(outunc, 3)
    nexttile

    % Add zero line to subplots 2, 3, and 4
    if i > 1
        plot(TT, zeros(size(TT)), 'k--', 'LineWidth', 0.5);
        hold on;
    end
    
    % Plot the fan chart for 'outunc'
    h_outunc = plot(TT, median(outunc(:, end-size(TT,2)+1:end, i)), 'Color',fan_color, 'LineWidth', 1.5);
    hold on;
    fill([TT, fliplr(TT)], [pct_values_outunc{i}(1,:), fliplr(pct_values_outunc{i}(end,:))],fan_color, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    for j = 2:numel(percentiles)-1
        fill([TT, fliplr(TT)], [pct_values_outunc{i}(j,:), fliplr(pct_values_outunc{i}(end-j+1,:))],fan_color, 'EdgeColor', 'none', 'FaceAlpha', 0.075);
    end
    axis tight;

    % Add shaded grey area for the period of forecast (sent to the background)
    forecast_start = datetime(2023, 03, 01);
    forecast_end = datetime(2024, 03, 01);
    forecast_area = [forecast_start forecast_end forecast_end forecast_start];
    y_limits = ylim;
    h_area = fill(forecast_area, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    uistack(h_area,'bottom')

    xticks(TT(1:3:end));
    xtickformat('MMM-yy')
    
    hold off; 
    set(gca, 'Layer', 'Top','FontSize',8);
    title(titles{i},"FontSize",12);
end

% Add legend
leg = legend( h_outunc, {'Unconditional Forecast (Median)'},'Orientation', 'horizontal','FontSize', 12);
leg.Layout.Tile = 'south';

set(gcf,'position',[x0,y0,width,height])
exportgraphics(t, 'unconditional_forecast.png', 'BackgroundColor','white','Resolution', 300);


