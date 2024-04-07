% function Plot_gSM(KM1,KM2,SMo)
function Plot_gSM(FC,WP,SMo)

close all
% Define the function
% g_SM_calc = @(KM1, KM2, SM, SMo) max(1 - KM1 * exp(KM2 * (SM - SMo)), 0);
g_SM_calc = @(FC, WP, SM)   min( max( (SM-WP)./(FC-WP), 0), 1);

% Define the range of SM/SMo values
SM_SMo_ratio = linspace(0, 2, 100);  % Adjust the range as needed


% Calculate g_SM values
SM   = 1:1:SMo;
% g_SM = g_SM_calc(KM1, KM2, SM, SMo);
g_SM = g_SM_calc(FC, WP, SM);

SM_SMo_ratio = SM./SMo;

% Plot
figure(1)
subplot(1,2,1)
plot(SM_SMo_ratio, g_SM, 'LineWidth', 2);
xlabel('SM / SMo');
ylabel('g_{SM}');
title('Plot of g_{SM} versus SM / SMo');
grid on;
% KM1 = 3.36*10^-4;  % Example values, adjust as needed
axis([ 0 1 0 1])
subplot(1,2,2)
plot(SM, g_SM, 'LineWidth', 2);
xlabel('SM');
ylabel('g_{SM}');
grid on;
