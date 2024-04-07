clear all
clc 
close all

%% 1 - Load Data
Q              = csvread('Q_mm.csv');   % Q in mm 
Dates          = Q(:,1:3); 
Precip         = csvread('Precip.csv'); % P in mm 
e              = csvread('e.csv');      % e in kPa
u2             = csvread('u2.csv');     % Wind Speed at 2m in m/s
S_in           = csvread('S_in.csv');   % Solar radiation in W/m2
Temp           = csvread('Temp.csv');   % Temp in Celsius
Lat_Lon_A_Z    = csvread('Lat_Lon_Area_Z.csv');   % Elevation

Q(:,1:3)       = [];
Precip(:,1:3)  = [];
e(:,1:3)       = [];
u2(:,1:3)      = [];
S_in(:,1:3)    = [];
Temp(:,1:3)    = [];


%% 2.1 - Generate Daily PET: Calculate Net Longwave

[Ra,Rso]    = Rso_calc(Dates,Lat_Lon_A_Z(:,1),Lat_Lon_A_Z(:,2),Lat_Lon_A_Z(:,4));    % Rso = Clear sky solar radiation MJ/(m2.day);

Rso        = Rso./(10^-6*60*60*24); % Convert to W/m2
Ra         = Ra./(10^-6*60*60*24); % Convert to W/m2

from  = 100; % start plot at
to    = 300; % end plot at
chose = 15;  % pick catchment

figure(1)
clf(1)
plot(Ra(from:to,chose),'.b'); hold on
plot(Rso(from:to,chose),'.k'); hold on
plot(S_in(from:to,chose),'.-r'); hold on
legend('Extraterrestrial solar radiation','Clear sky solar radiation','S in')
ylabel('W/m^2')
xlabel('days')

f        = S_in./Rso;
f        = 1.35.*(S_in./Rso) - 0.35;

e_prime  = 0.34 - 0.14*(e./1000).^0.5;
L_net    = -f.*e_prime.*5.67*10^-8.*(Temp+273.15).^4; % Equation 5.22 Terrestrial Hydrometeorology (TH)


%% 2.2 - Generate Daily PET:
clc
Penman_Monteith_calc = @(Delta, A, rho_a, c_p, VPD, ra, gamma, rs)          ( Delta.*A + rho_a.*c_p.*VPD./ra ) ./ ( Delta + gamma.*(1 + rs./ra) );


cp       = 1.1013.*10^3;                                  % specific heat at constant pressure for air J/(kg.K)
lambda   = 10^3*(2.501 - 0.002361.*Temp);                 % Latent Heat of Vaporization (kJ/kg)
rho_air  = 1.23;                                          % air density kg/m3

ra       = 208./u2;
rs       = 70; 
Rn_RC    = S_in.*(1-0.23) + L_net;                        % Net Radiation over a Reference Crop, in W/m2
e_sat    = 0.6108.*exp( (17.27.*Temp)./(237.3 + Temp) );  % Saturated Vapor Pressure in kPa
D        = e_sat - e;                                     % Vapor Pressure Deficit in kPa
Delta    = 4098.*e_sat./(237.3 + Temp).^2;                % Slope of esat versus temp curve (kPa/C)
Z        = Lat_Lon_A_Z(:,4);                              % Elevation (m)
Press    = 101.3*((293-0.0065.*Z)/293).^5.26;             % Pressure as function of elevation (kPa)
gamma    = cp.*(Press)'./(0.622.*lambda.*10^3);           % Psychrometric Constant (kPa/degC)

Wm2_mm  = (lambda.*10^3).^-1.*86400;                      % Conversion Factor from W/m2 to mm/day


E_RC           = ( Delta.*Rn_RC + rho_air.*cp.*D./ra )./( Delta + gamma.*( 1 + rs./ra) );
E_RC_backup    = Penman_Monteith_calc(Delta,Rn_RC,rho_air,cp,D,ra,gamma,rs);

E_RC_mm  = E_RC.*Wm2_mm;
Rn_RC_mm = Rn_RC.*Wm2_mm;

gamma_m   = gamma.*(1 + 0.33.*u2) ;
% E_RC_mm_2 = ( Delta./(Delta + gamma_m) ).*Rn_RC.*Wm2_mm + ( gamma_m./(Delta + gamma_m) ).*( 900./(Temp + 273) ).*u2.*D;
%% 2.3 Pick a catchment

% pick_catchment = 2;
% pick_catchment = 9;
% pick_catchment = 8;
pick_catchment = 10;

%% 3.1 Mean annual analysis
n_years = length(unique(Dates(:,1)));

[PET_MEAN_MONTHLY,PET_MEAN_ANNUAL]       = make_means_new(E_RC_mm,Dates);
[Q_MEAN_MONTHLY,Q_MEAN_ANNUAL]           = make_means_new(Q,Dates);
[P_MEAN_MONTHLY,P_MEAN_ANNUAL]           = make_means_new(Precip,Dates);
[Temp_MEAN_MONTHLY,Temp_MEAN_ANNUAL]     = make_means(Temp,Dates);


PHI             = PET_MEAN_ANNUAL./P_MEAN_ANNUAL;
E_P             = (P_MEAN_ANNUAL-Q_MEAN_ANNUAL)./P_MEAN_ANNUAL;
% Budyko plot
AI_plot_lin         = 0.000:0.01:50;
Budyko              = @(AI)    (AI.*(1-exp(-AI)).*tanh(AI.^-1)).^0.5;

figure(5)
clf(5)
set(gcf,"Position", [    488.0000  548.2000  724.2000  213.8000])
subplot(1,2,1)
plot(PHI,E_P, 'ok', 'LineWidth',1);hold on
plot(PHI(pick_catchment),E_P(pick_catchment), '.r', 'MarkerSize',20);hold on

plot(AI_plot_lin,Budyko(AI_plot_lin),'k-','Linewidth',1);hold on

legend(' Catchments',' Your pick' ,' Budyko')
legend boxoff; ylabel('E/P'); xlabel('Phi')
axis([ 0 5 0 1])

subplot(1,2,2)
plot([1:12],P_MEAN_MONTHLY(:,pick_catchment), '-b', 'LineWidth',1); hold on
plot([1:12],PET_MEAN_MONTHLY(:,pick_catchment), '-r', 'LineWidth',1); hold on
ylabel('mm')

yyaxis right
plot([1:12],Temp_MEAN_MONTHLY(:,pick_catchment), '--m', 'LineWidth',1); hold on
ylabel('deg. C')
xlabel('Month')
legend('P','PET', 'Temp')
%% 3.2 Run Hydrologic Model
clc
clear INPUT;

INPUT(:,1)  = Precip(:,pick_catchment);
INPUT(:,2)  = E_RC_mm(:,pick_catchment);
INPUT(:,3)  = e(:,pick_catchment);
INPUT(:,4)  = u2(:,pick_catchment);
INPUT(:,5)  = S_in(:,pick_catchment);
INPUT(:,6)  = Temp(:,pick_catchment);
INPUT(:,7)  = L_net(:,pick_catchment);
Dates_sim   =  Dates(:,:);

vars_forest   =  water_energy_variables(1);
vars_grass    =  water_energy_variables(2);


%% 3.3 Simulation 1
PAR(1) = 50;     % % Maximum infiltration rate  [300]
PAR(2) = 150;   % Su_max - total water capacity [50-300]
PAR(3) = 30;    % Ts time - parameter for slowflow [20-100]
PAR(4) = 1;     % Tf time = parameter for quickflow [1-3]

vars_1            = vars_forest;
vars_1.Patm       = Press(pick_catchment);
vars_1.g0         = 10;
vars_1.z0         = 0.826;
vars_1.a          = 0.12;

Su_max                   = PAR(2);
vars_1.FC                = 0.35.*Su_max;
vars_1.WP                = 0.11.*Su_max;

Plot_gSM(vars_1.FC,vars_1.WP,Su_max)

%
OUT_1             = toymodel_update(INPUT, PAR, vars_1);
OUT_1             = aggregate(OUT_1, Dates_sim);


% Check your results
PHI_SIM            = OUT_1.PET_mean_annual./OUT_1.P_mean_annual;
EF_SIM             = (OUT_1.P_mean_annual - OUT_1.QT_mean_annual)./OUT_1.P_mean_annual;
% EF_SIM             = OUT_1.E_mean_annual./OUT_1.P_mean_annual;

Q_SIM_mean_monthly = OUT_1.QT_mean_mon;
Q_SIM              = OUT_1.QT;
plot_from          = 1000;
plot_to            = 1300;

figure(5)
clf(5)
set(gcf,"Position", 1.0e+03 *[0.2490    0.5210    1.0696    0.2410])

subplot(1,2,1)
plot(PHI,E_P, 'ok', 'LineWidth',1);hold on
plot(PHI(pick_catchment),E_P(pick_catchment), '.r', 'MarkerSize',20);hold on
plot(PHI_SIM,EF_SIM, '.b', 'MarkerSize',20);hold on
plot(AI_plot_lin,Budyko(AI_plot_lin),'k-','Linewidth',1);hold on
legend(' Catchments',' Your catchment' ,' Your catchment (simulated)',' Budyko' ,'Location', 'east'); 
legend boxoff; ylabel('E/P'); xlabel('Phi')
axis([ 0 5 0 1])

subplot(1,2,2)
plot([1:12],Q_SIM_mean_monthly, '--k', 'LineWidth',1);hold on
plot([1:12],Q_MEAN_MONTHLY(:,pick_catchment), '--r', 'LineWidth',1);hold on
legend(' Simulated', 'Observed'); 


figure(6)
clf(6)
set(gcf,"Position", [320.2000  105.8000  830.4000  312.0000])
subplot(1,2,1)

plot(Q_SIM(plot_from:plot_to), '-k', 'LineWidth',0.5);hold on
plot(Q(plot_from:plot_to,pick_catchment), '-r', 'LineWidth',0.5);hold on
legend(' Simulated', 'Observed'); 

subplot(1,2,2)
plot(OUT_1.Su, '-k', 'LineWidth',0.5);hold on
legend(' Simulated');
yyaxis right
% plot(OUT_1.R, '-r', 'LineWidth',0.5);hold on
plot(OUT_1.QS, '-r', 'LineWidth',0.5);hold on

%% 3.4 Simulation 2
vars_2A            = vars_forest;
vars_2A.Patm       = Press(pick_catchment);
vars_2A.g0         = 10;
vars_2A.z0         = 0.83;
vars_2A.a          = 0.12;
vars_2A.FC          = 0.38.*Su_max;
vars_2A.WP          = 0.11.*Su_max;


vars_2B            = vars_grass;
vars_2B.Patm       = Press(pick_catchment);
vars_2B.g0         = 3;
vars_2B.z0         = 0.24;
vars_2B.a          = 0.25;
vars_2B.FC          = 0.38.*Su_max;
vars_2B.WP          = 0.11.*Su_max;
PAR2B               = PAR;
PAR2B(1)            = 25;

weight_A           = 0.3;
weight_B           = 0.7;

OUT_2A             = toymodel_update(INPUT, PAR, vars_2A);
OUT_2B             = toymodel_update(INPUT, PAR2B, vars_2B);

OUT_2              = weighted_average_outputs(OUT_2A, OUT_2B, weight_A, weight_B);
OUT_2              = aggregate(OUT_2, Dates_sim);



%% 3.5 Produce final figures:



