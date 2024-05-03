clear all
clc
close all
%% Part 0: Load Data
Lower_Hafren = csvread("data.csv",1,0);
P    = Lower_Hafren(:,3);
Q    = Lower_Hafren(:,5);
PET  = Lower_Hafren(:,6)/0.7;
Cl_P = Lower_Hafren(:,4);

excelSerialDates       = Lower_Hafren(:,2);
mdates                 = excelSerialDates + datenum('1899-12-30');

Dates(:,1)             = year(mdates);
Dates(:,2)             = month(mdates);
Dates(:,3)             = day(mdates);
%% Part 0: Load Data
Cl             = Lower_Hafren(:,7);
Cl(Cl<0.00001) = NaN;
% Find indices of NaN values
nanIndices = isnan(Cl);

% Find indices of non-NaN values
nonNanIndices = find(~isnan(Cl));

% Create index vector for NaN values
nanIndexVector = 1:numel(Cl);
nanIndexVector = nanIndexVector(nanIndices);

% Interpolate NaN values based on non-NaN values
ClInterpolated = interp1(nonNanIndices, Cl(nonNanIndices), nanIndexVector, 'linear');

% Replace NaN values with interpolated values
Cl(nanIndexVector) = ClInterpolated;

%% Part 0: Observed Data
plot_from = find(Dates(:,1) == 1988,1,'first');
plot_to   = find(Dates(:,1) == 1990,1,'last');

figure(1)
clf(1)
subplot(1,2,1)
a = plot(mdates(plot_from:plot_to),P(plot_from:plot_to),'k-', LineWidth=0.5); hold on
b = plot(mdates(plot_from:plot_to),Q(plot_from:plot_to),'-b', LineWidth=0.5); hold on
legend([a,b],'P','Q');ylabel('mm/day')
datetick('x', 'mm/dd/yyyy', 'keeplimits')

subplot(1,2,2)
c = plot(mdates(plot_from:plot_to),Cl_P(plot_from:plot_to),'r-'); hold on
d = plot(mdates(plot_from:plot_to),Cl(plot_from:plot_to),'ko', Markersize=3); hold on

set(gca, 'YColor', 'k');ylim([ 0 30]);
legend([c,d],'Cl_P','Cl_Q');ylabel('Cl (ppm)')
datetick('x', 'mm/dd/yyyy', 'keeplimits')

%% Part 0: The simplest hydro model possible ( A linear reservoir)

Su_max    = 200;
K         = 3;
Su        = zeros(length(Dates),1);
ET        = zeros(length(Dates),1);
Q_sim     = zeros(length(Dates),1);
Su_begin  = 1000;
Su_old    = Su_begin;

for i= 1:length(Dates)

ET(i,1)    = PET(i,1)*(Su_old/Su_max);
Q_sim(i,1) = Su_old*(1/K);

Su(i,1)    = Su_old + P(i,1) - ET(i,1) - Q_sim(i,1); 
Su_old     = Su(i,1);

end

figure(2)
clf(2)
a = plot(mdates(plot_from:plot_to), Q(plot_from:plot_to), '-', 'Color', [0.3, 0.3, 0.3]); 
hold on
b = plot(mdates(plot_from:plot_to), Q_sim(plot_from:plot_to), 'b-'); 
hold on
ylabel('Q (mm)')
yyaxis right
c = plot(mdates(plot_from:plot_to), ET(plot_from:plot_to), 'Color', [0, 0.6, 0]); 
hold on
ylabel('ET (mm)')
set(gca, 'YColor', [0, 0.6, 0]); % Set right side y-axis color to match plot c
ylim([-5 5]); % Uncomment this line if you want to set custom y-axis limits
legend([a, b, c], 'Q_{obs}', 'Q_{sim}', 'ET_{sim}')
datetick('x', 'mm/yyyy', 'keeplimits')

%% Part 1: Define the initial age structure
mu               = 250;                   % Mean age
sigma            = 200;                   % Standard deviation
vector_length    = length(Q);             % Length of the vector
S_init           = 400;                   % Initial Storage in mm
ages_initial_length = 1000;
% dummy            = normpdf([1:1:ages_initial_length], mu, sigma);
% dummy            = exppdf([1:1:ages_initial_length], mu);
dummy              = ones(1, ages_initial_length) / ages_initial_length;

T_pdf_temp       = NaN(vector_length, 1);
T_pdf_temp(1:ages_initial_length) = dummy;
s_T_temp         = T_pdf_temp.*S_init;
S_T_temp         = cumsum(s_T_temp);

S_T(:, 1)        = S_T_temp; 
s_T(:, 1)        = s_T_temp; 

figure(3)
clf(3)
subplot(1,3,1)
plot(T_pdf_temp);title('Chosen pdf'); xlabel('Water Ages'); ylabel('Prob. (1/age)');
subplot(1,3,2)
plot(s_T_temp);title('initial s_T'); xlabel('Water Ages'); ylabel('Storage. (mm/age)');
subplot(1,3,3)
plot(S_T_temp);title('initial S_T'); xlabel('Water Ages'); ylabel('Cumulative Storage. (mm/age)');
%% Part 1: Define Simulation length

sim_length = 3000;
%% Part 2: Water ages over time
clear S_T  s_T S_T_new s_T_new S_T  

S_T(:, 1)        = S_T_temp; 
s_T(:, 1)        = s_T_temp; 
for i = 2:sim_length

S_T_new      = S_T(:,i-1);


% Ageing and adding Water with age=0 
s_T_new        = [S_T_new(1); diff(S_T_new)];
s_T_new        = [0; s_T_new];
s_T_new(end)   = [];

S_T_new        = cumsum(s_T_new);

S_T(:,i)         = S_T_new;
s_T(:,i)         = s_T_new;



end
STORAGE = nansum(s_T);

%Plot aging
clear time_stamps_ch
time_stamps        = [1,50,200,500,1000];

for i =1:5
time_stamps_ch{i} = num2str(time_stamps(i));
end

colors = cool(length(time_stamps)); % Choosing colors from 'winter' colormap

figure(4)
clf(4)
set(gcf,"Position",1.0e+03 *[ 0.3802    0.4906    1.0416    0.271])
for i = 1:length(time_stamps)
    subplot(1,3,1)
    plot(s_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('s_T')

    subplot(1,3,2)
    plot(S_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('S_T')
end
legend(time_stamps_ch); hold on
subplot(1,3,3)
plot(STORAGE, 'Color', 'black'); hold on; title('Storage (mm)');
ylim([0 S_init*2])
%% Part 3: Water ages over time + new water is added.
clear S_T  s_T S_T_new s_T_new S_T  

S_T(:, 1)        = S_T_temp; 
s_T(:, 1)        = s_T_temp; 
s_T_retrieve     = [S_T(1); diff(S_T_temp)];

for i = 2:sim_length
                                       
   

    S_T_new      = S_T(:,i-1);


% Ageing and adding Water with age=0 
s_T_new        = [S_T_new(1); diff(S_T_new)];
s_T_new        = [0; s_T_new];
s_T_new(1)     = P(i);
s_T_new(end)   = [];

S_T_new            = cumsum(s_T_new);

S_T(:,i)         = S_T_new;
s_T(:,i)         = s_T_new;

end
STORAGE = nansum(s_T);


%Plot aging
clear time_stamps_ch
time_stamps        = [1,50,200,500,1000];

for i =1:5
time_stamps_ch{i} = num2str(time_stamps(i));
end

colors = cool(length(time_stamps)); % Choosing colors from 'winter' colormap

figure(5)
clf(5)
set(gcf,"Position",1.0e+03 *[ 0.3802    0.4906    1.0416    0.271])
for i = 1:length(time_stamps)
    subplot(1,3,1)
    plot(s_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('s_T')

    subplot(1,3,2)
    plot(S_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('S_T')
end
legend(time_stamps_ch); hold on
subplot(1,3,3)
plot(STORAGE, 'Color', 'black'); hold on; title('Storage (mm)');


%% Part 4: OMEGA Q as a gamma 
% close all
% mean_OMEGA_Q = 10;
% var_OMEGA_Q  = 10;
% 
% k_gamma     = (mean_OMEGA_Q^2) / var_OMEGA_Q;
% theta_gamma = var_OMEGA_Q / mean_OMEGA_Q;
% 
% 
% % Calculate the theoretical PDF and CDF of gamma distribution
% x_gamma     = linspace(0, ages_initial_length, ages_initial_length);
% omega_Q     = gampdf(x_gamma, k_gamma, theta_gamma);
% OMEGA_Q     = gamcdf(x_gamma, k_gamma, theta_gamma);
% 
% figure(1)
% clf(1)
% subplot(1,2,1)
% plot(x_gamma, omega_Q)
% title('\omega_{Q} (Gamma Distribution)')
% subplot(1,2,2)
% plot(x_gamma, OMEGA_Q)
% title('\Omega_{Q} (Gamma Distribution)')

%% Part 4: OMEGA Q as exponential
mean_OMEGA_Q    = 100;
mean_exp        = mean_OMEGA_Q;

% Calculate the theoretical PDF and CDF of exponential distribution
x_exp     = linspace(0, ages_initial_length, ages_initial_length);
omega_Q     = exppdf(x_exp, mean_exp);
OMEGA_Q     = expcdf(x_exp, mean_exp);

figure(1)
clf(1)
subplot(1,2,1)
plot(x_exp, omega_Q)
title('\omega_{Q} (Gamma Distribution)')
subplot(1,2,2)
plot(x_exp, OMEGA_Q)
title('\Omega_{Q} (Gamma Distribution)')


%% Part 4: Water ages over time + new water is added + removal through Q
clear S_T  s_T S_T_new s_T_new S_T   OMEGA_Q_SAVE Su_actual Su_only_P

S_T(:, 1)        = S_T_temp; 
s_T(:, 1)        = s_T_temp; 
s_T_retrieve     = [S_T(1); diff(S_T_temp)];

Q_sim_actual        = Q_sim(28:end);
P_actual            = P(28:end);
Su_begin            = nansum(s_T);
Su_actual(1,1)      = Su_begin;
Su_only_P(1,1)      = Su_begin;

for i = 2:sim_length
                                       
                                        %Estimate the volume of water leaving as streamflow with varying ages (in CDF format)
                                        x_OMEGA     = linspace(0, ages_initial_length, ages_initial_length+i);
%                                         OMEGA_Q     = gamcdf(x_OMEGA, k_gamma, theta_gamma);
                                        OMEGA_Q     = expcdf(x_OMEGA, mean_exp);

                                        out_q_aux                                =  Q_sim_actual(i).*OMEGA_Q';
                                        out_q                                    =  NaN(vector_length,1);
                                        out_q(1:length(OMEGA_Q))                 =  out_q_aux;
                                        out_q(length(OMEGA_Q):end)               =  out_q_aux(length(OMEGA_Q));
                
                                          % Remove storage with varying ages from S_T using OMEGA functions          
                                         [S_T_new, out_q]    = correct_S_T_new(S_T(:,i-1),out_q); OMEGA_Q_SAVE(:,i)  =  out_q;


s_T_new              = [S_T_new(1); diff(S_T_new)];

% Ageing and adding Water with age = 0 
s_T_new        = [0; s_T_new];
s_T_new(1)     = P_actual(i);
s_T_new(end)   = [];

S_T_new            = cumsum(s_T_new);

S_T(:,i)         = S_T_new;
s_T(:,i)         = s_T_new;



% compute storage variations
Su_actual(i,1)         = Su_actual(i-1,1) + P_actual(i) - Q_sim_actual(i);
Su_only_P(i,1)         = Su_only_P(i-1,1) + P_actual(i);

end
STORAGE = nansum(s_T);



%% Part 4: Plots
clear time_stamps_ch
time_stamps        = [1,1000,1500,2000,3000];

for i =1:5
time_stamps_ch{i} = num2str(time_stamps(i));
end

colors    = cool(length(time_stamps)); % Choosing colors from 'winter' colormap
xplot_min = 0;
xplot_max = 1000;

figure(6)
clf(6)
set(gcf,"Position",1.0e+03 *[ 0.3802    0.4906    1.0416    0.271])
for i = 1:length(time_stamps)
    subplot(1,4,1)
    plot(s_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('s_T')
    xlim([xplot_min xplot_max])

    subplot(1,4,2)
    plot(S_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('S_T')
    xlim([xplot_min xplot_max])

    subplot(1,4,3)
    plot(OMEGA_Q_SAVE(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('\Omega_Q x Q(t)')
    xlim([xplot_min xplot_max])

end
legend(time_stamps_ch); hold on
subplot(1,4,4)
plot(STORAGE, 'Color', 'black'); hold on; title('Storage (mm)');

%% Part 5: OMEGA ET as uniform

% Define parameters for the uniform distribution
a_ET = 0;  % Lower bound
b_ET = 10; % Upper bound

% Generate data for omega_ET
x_uniform = linspace(a_ET, b_ET, b_ET-a_ET);
omega_ET  = unifpdf(x_uniform, a_ET, b_ET);

% Calculate cumulative distribution function of omega_ET
OMEGA_ET = unifcdf(x_uniform, a_ET, b_ET);


figure(1)
clf(1)

subplot(1,2,1)
plot(x_uniform, omega_ET)
title('\omega_{ET} (Uniform Distribution)')
subplot(1,2,2)
plot(x_uniform, OMEGA_ET)
title('\Omega_{ET} (Uniform Distribution)')
%% Part 5: Water ages over time + new water is added + removal through Q + removal with ET
clear S_T  s_T S_T_new s_T_new S_T   Su_only_P Su_only_P OMEGA_Q_SAVE OMEGA_ET_SAVE  Cl_S

S_T(:, 1)           = S_T_temp; 
s_T(:, 1)           = s_T_temp; 
Cl_S(:, 1)          = NaN(vector_length,1); 
Cl_S_old            = Cl_S(:, 1);

Q_sim_actual        = Q_sim(28:end);
% Q_sim_actual        = Q(28:end);

P_actual            = P(28:end);
E_actual            = ET(28:end);
Cl_P_actual         = Cl_P(28:end);

Su_begin            = nansum(s_T);
Su_actual(1,1)      = Su_begin;
Su_only_P(1,1)      = Su_begin;

for i = 2:sim_length
                                       
                                        %Estimate the volume of water leaving as streamflow with varying ages (in CDF format)
                                        x_OMEGA     = linspace(0, ages_initial_length, ages_initial_length+i);
%                                         OMEGA_Q     = gamcdf(x_OMEGA, k_gamma, theta_gamma);
                                        OMEGA_Q     = expcdf(x_OMEGA, mean_exp);

                                        out_q_aux                                =  Q_sim_actual(i).*OMEGA_Q';
                                        out_q                                    =  NaN(vector_length,1);
                                        out_q(1:length(OMEGA_Q))                 =  out_q_aux;
                                        out_q(length(OMEGA_Q):end)               =  out_q_aux(length(OMEGA_Q));
                
                                        %Estimate the volume of water leaving as ET with varying ages (in CDF format)
                                        OMEGA_ET                                 = unifcdf(x_uniform, a_ET, b_ET);      % ASSUMES ET ALWAYS SAMPLE FROM A FIXED AGE POOL
                                        out_et_aux                               =  E_actual(i).*OMEGA_ET';
                                        out_et                                   =  NaN(vector_length,1);
                                        out_et(1:length(OMEGA_ET))               =  out_et_aux;
                                        out_et(length(OMEGA_ET):end)             =  out_et_aux(length(OMEGA_ET));


                                          % Remove storage with varying ages from S_T using OMEGA functions 
                                         [S_T_new, out_et]    = correct_S_T_new(S_T(:,i-1),out_et); OMEGA_ET_SAVE(:,i)  =  out_et;
                                         [S_T_new, out_q]     = correct_S_T_new(S_T_new,out_q);     OMEGA_Q_SAVE(:,i)   =  out_q;


s_T_new              = [S_T_new(1); diff(S_T_new)];
% Ageing and adding Water with age = 0 
s_T_new        = [0; s_T_new];
s_T_new(1)     = P_actual(i);
s_T_new(end)   = [];

S_T_new            = cumsum(s_T_new);

S_T(:,i)         = S_T_new;
s_T(:,i)         = s_T_new;

% Keep track of Chloride concentrations added to storage
Cl_S_new         = [Cl_P_actual(i); Cl_S_old];
Cl_S_new(end)    = [];
Cl_S(:,i)        = Cl_S_new;     
Cl_S_old         = Cl_S_new;

% compute storage variations
Su_actual(i,1)         = Su_actual(i-1,1) + P_actual(i) - Q_sim_actual(i);
Su_only_P(i,1)         = Su_only_P(i-1,1) + P_actual(i);

end
STORAGE = nansum(s_T);

%% Part 5: Plots
clear time_stamps_ch
time_stamps        = [1,1000,1500,2000,3000];

for i =1:5
time_stamps_ch{i} = num2str(time_stamps(i));
end

colors    = cool(length(time_stamps)); % Choosing colors from 'winter' colormap
xplot_min = 0;
xplot_max = 1000;

figure(6)
clf(6)
set(gcf,"Position",1.0e+03 *[ 0.3802    0.4906    1.0416    0.271])
for i = 1:length(time_stamps)
    subplot(1,5,1)
    plot(s_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('s_T')
    xlim([xplot_min xplot_max])

    subplot(1,5,2)
    plot(S_T(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('S_T')
    xlim([xplot_min xplot_max])

    subplot(1,5,3)
    plot(OMEGA_Q_SAVE(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('\Omega_Q x Q(t)')
    xlim([xplot_min xplot_max])

    subplot(1,5,4)
    plot(OMEGA_ET_SAVE(:,time_stamps(i)), 'Color', colors(i, :)); hold on; title('\Omega_{ET} x ET(t)')
    xlim([xplot_min xplot_max])
end
legend(time_stamps_ch); hold on
subplot(1,5,5)
plot(STORAGE, 'Color', 'black'); hold on; title('Storage (mm)');


%% Part 6: Tracer concentrations at the outflows:

clear mass_cons Cl_Q Cl_ET

for i = 3000:-1:2
             
    omega_Q       = [OMEGA_Q_SAVE(i,1);  diff( OMEGA_Q_SAVE(:,i) )  ] ;
    omega_ET      = [OMEGA_ET_SAVE(i,1);  diff( OMEGA_ET_SAVE(:,i) )  ] ;

    Cl_history   = Cl_S(:,i-1);

    %%%% Make sure mass is being converved
    vec_aux           = Cl_history./Cl_history;
    test              = omega_Q.*vec_aux;
    mass_cons(i,1)    = nansum(test)-nansum(omega_Q);
    %%%% 

    
    Cl_Q(i,1)          = nansum( omega_Q./nansum(omega_Q).*Cl_history );
    Cl_ET(i,1)         = nansum( omega_ET./nansum(omega_ET).*Cl_history );

end

% Part 6: Plot
figure(8)
clf(8)
set(gcf,'Position',1.0e+03*[0.3098    0.4018    1.0584    0.3602])
mdates_range = plot_from:plot_to;
a = plot(mdates(mdates_range), Cl_P(mdates_range), 'r-', 'LineWidth', 0.5); 
hold on
b = plot(mdates(mdates_range), Cl(mdates_range), 'ko', 'MarkerSize', 3); 
hold on
c = plot(mdates(mdates_range), Cl_Q(mdates_range - 28), 'bo', 'MarkerSize', 3); 
hold on
d = plot(mdates(mdates_range), Cl_ET(mdates_range - 28), 'o', 'Color', [0, 0.6, 0], 'MarkerSize', 3); 
hold on
set(gca, 'YColor', 'k');
ylim([0 40]);
legend([a,b,c, d], 'Cl_P', 'Cl_Q (observed)',' Cl_Q (simulated)',' Cl_ET (simulated)');
ylabel('Cl (ppm)');
datetick('x', 'mm/dd/yyyy', 'keeplimits');

% Set the range of the x-axis
xlim([mdates(plot_from) mdates(plot_to)]);
%% Part 7: OMEGA Q as exponential with storage-dependent mean
close all
% Define the range of x values
x = linspace(0, 1, 100); % 100 points from 0 to 1

% Define the maximum and minimum values of Y
max_mean_exp = 10; % Maximum value of Y
min_mean_exp = 1;  % Minimum value of Y

% Calculate Y using a linear interpolation between max_Y and min_Y
Y = max_mean_exp - (max_mean_exp - min_mean_exp) * x;

Su_ratio          = Su./Su_max;
mean_exp_variable =  max_mean_exp - (max_mean_exp - min_mean_exp) * Su_ratio;

mean_exp_variable(mean_exp_variable<min_mean_exp)=min_mean_exp;

% Plotting
figure(1)
clf(1)
subplot(1,2,1)
plot(x, Y, 'LineWidth', 1);
xlabel('Su/Sumax');
ylabel('mean_{exp}');
title('Variation of Y with x');
grid on;
subplot(1,2,2)
plot(Su(mdates_range), 'LineWidth', 1);
yyaxis right
plot(mean_exp_variable(mdates_range), 'LineWidth', 1);
legend('Su','Mean of exp.')


%% Part 7: Time varying selection


clear S_T  s_T S_T_new s_T_new S_T   Su_only_P Su_only_P OMEGA_Q_SAVE OMEGA_ET_SAVE  Cl_S

S_T(:, 1)           = S_T_temp; 
s_T(:, 1)           = s_T_temp; 
Cl_S(:, 1)          = NaN(vector_length,1); 
Cl_S_old            = Cl_S(:, 1);

Q_sim_actual        = Q_sim(28:end);
Q_sim_actual        = Q(28:end);

P_actual            = P(28:end);
E_actual            = ET(28:end);
Cl_P_actual         = Cl_P(28:end);
mean_exp_var_actual = mean_exp_variable(28:end);

Su_begin            = nansum(s_T);
Su_actual(1,1)      = Su_begin;
Su_only_P(1,1)      = Su_begin;

for i = 2:sim_length
                                       
                                        %Estimate the volume of water leaving as streamflow with varying ages (in CDF format)
                                        x_OMEGA     = linspace(0, ages_initial_length, ages_initial_length+i);
%                                         OMEGA_Q     = gamcdf(x_OMEGA, k_gamma, theta_gamma);
                                        OMEGA_Q     = expcdf(x_OMEGA, mean_exp_var_actual(i));

                                        out_q_aux                                =  Q_sim_actual(i).*OMEGA_Q';
                                        out_q                                    =  NaN(vector_length,1);
                                        out_q(1:length(OMEGA_Q))                 =  out_q_aux;
                                        out_q(length(OMEGA_Q):end)               =  out_q_aux(length(OMEGA_Q));
                
                                        %Estimate the volume of water leaving as ET with varying ages (in CDF format)
                                        OMEGA_ET                                 = unifcdf(x_uniform, a_ET, b_ET);      % ASSUMES ET ALWAYS SAMPLE FROM A FIXED AGE POOL
                                        out_et_aux                               =  E_actual(i).*OMEGA_ET';
                                        out_et                                   =  NaN(vector_length,1);
                                        out_et(1:length(OMEGA_ET))               =  out_et_aux;
                                        out_et(length(OMEGA_ET):end)             =  out_et_aux(length(OMEGA_ET));


                                          % Remove storage with varying ages from S_T using OMEGA functions 
                                         [S_T_new, out_et]    = correct_S_T_new(S_T(:,i-1),out_et); OMEGA_ET_SAVE(:,i)  =  out_et;
                                         [S_T_new, out_q]     = correct_S_T_new(S_T_new,out_q);     OMEGA_Q_SAVE(:,i)   =  out_q;


s_T_new              = [S_T_new(1); diff(S_T_new)];
% Ageing and adding Water with age = 0 
s_T_new        = [0; s_T_new];
s_T_new(1)     = P_actual(i);
s_T_new(end)   = [];

S_T_new            = cumsum(s_T_new);

S_T(:,i)         = S_T_new;
s_T(:,i)         = s_T_new;

% Keep track of Chloride concentrations added to storage
Cl_S_new         = [Cl_P_actual(i); Cl_S_old];
Cl_S_new(end)    = [];
Cl_S(:,i)        = Cl_S_new;     
Cl_S_old         = Cl_S_new;

% compute storage variations
Su_actual(i,1)         = Su_actual(i-1,1) + P_actual(i) - Q_sim_actual(i);
Su_only_P(i,1)         = Su_only_P(i-1,1) + P_actual(i);

end
STORAGE = nansum(s_T);
%% Part 7: Plot

[Cl_Q_2,Cl_ET_2] = get_conc_out(Cl_S,OMEGA_Q_SAVE,OMEGA_ET_SAVE);

figure(9)
clf(9)
set(gcf,'Position',1.0e+03*[0.3098    0.4018    1.0584    0.3602])
mdates_range = plot_from:plot_to;
a = plot(mdates(mdates_range), Cl_P(mdates_range), 'r-', 'LineWidth', 0.5); 
hold on
b = plot(mdates(mdates_range), Cl(mdates_range), 'ko', 'MarkerSize', 3); 
hold on
c = plot(mdates(mdates_range), Cl_Q_2(mdates_range - 28), 'bo', 'MarkerSize', 3); 
hold on
d = plot(mdates(mdates_range), Cl_ET_2(mdates_range - 28), 'o', 'Color', [0, 0.6, 0], 'MarkerSize', 3); 
hold on
set(gca, 'YColor', 'k');
ylim([0 40]);
legend([a,b,c, d], 'Cl_P', 'Cl_Q (observed)',' Cl_Q (simulated)',' Cl_ET (simulated)');
ylabel('Cl (ppm)');
datetick('x', 'mm/dd/yyyy', 'keeplimits');

% Set the range of the x-axis
xlim([mdates(plot_from) mdates(plot_to)]);
