%% TOY MODEL

% function [Ea,QF,R,QS,QT,Sf,Su,Ss,St,AL,IE,SE,     Ei,Et,S_canopy,pot_inf,  ET_vars]= toymodel_update(DATA,PAR,vars)
function OUT = toymodel_update(DATA,PAR,vars)

%% Inputs
% DATA = INPUT_2;
[M,~]   = size(DATA);
P       = DATA(:,1); 
Ep      = DATA(:,2);

e       = DATA(:,3);
u2      = DATA(:,4);
Sin     = DATA(:,5);
Temp    = DATA(:,6);
Lnet    = DATA(:,7);

%% Parameters
    mir        = PAR(1);  % Maximum infiltration rate 
    Su_max     = PAR(2);     % Unsaturated zone storage capacity [50-300]
    Ts         = PAR(3);  % time parameter for slowflow [20-100]
    Tf         = PAR(4);  % time parameter for quickflow [1-3]
    beta       = 1;      % Split between Recharge and overland flow [SET IT TO = 1]


%% Initialization
QF     = zeros(M,1);
QS     = zeros(M,1);
SE     = zeros(M,1);

QT     = zeros(M,1);
R      = zeros(M,1);
Sf     = zeros(M,1);
Su     = zeros(M,1);
Ss     = zeros(M,1);
Ea     = zeros(M,1);
AL     = zeros(M,1);
IE     = zeros(M,1);
SE     = zeros(M,1);

Ei     = zeros(M,1);
Et     = zeros(M,1);
S_canopy     = zeros(M,1);
pot_inf     = zeros(M,1);

% Initial values
S0          = 200; % initial condition
Su_dt       = S0;
Ss_dt       = S0;
Sf_dt       = S0;
S_canopy_old= 0;

%% MAIN FLOW ROUTINE
    
for t=1:M
%     fprintf('time: %d\n', t);

             if P(t)>0
             AL(t) = 1 - (1 - exp(-P(t)/mir) )/(P(t)/mir);
             else
             AL(t) = 0;
             end
             
            % Unsaturated zone water balance    
            [r,se,Su_dt,   S_can,E_int_mm,E_tran_mm,E_total_mm,rain_pass,ra,rs] = SU_eq_mod(Su_dt,P(t),Su_max,AL(t),beta,  vars,  Sin(t),Lnet(t),Temp(t),u2(t),e(t),S_canopy_old);

            ET_vars.ra(t,1)     = ra;
            ET_vars.rs(t,1)     = rs;

            R(t)         = r;              % recharge
            Ea(t)        = E_total_mm;     % evaporation 
            Su(t)        = Su_dt;
            SE(t)        = se;             % Saturation Excess Overland Flow Input [NOT CONSIDERED WHEN beta = 1]

            Ei(t)        = E_int_mm;       % save interception
            Et(t)        = E_tran_mm;      % save transpiration
            S_canopy(t)  = S_can;
            S_canopy_old = S_can;
            pot_inf(t)   = rain_pass;
            % Saturated zone
            [Qs,Ss_dt] = SS_eq(Ss_dt,r,Ts);
            
            QS(t) = Qs;    % slowflow
            Ss(t) = Ss_dt;  
end


for t=1:M
    
        
    [Qf,Sf_dt,ie]   = SF_eq(Sf_dt,P(t),AL(t),Tf,SE(t)); 

    QF(t) = Qf;         % quickflow
    Sf(t) = Sf_dt;      
    IE(t) = ie;

end

% model's outputs

QT   = QF + QS;
St   = Su + Ss;

 % Store outputs in the structure OUT
OUT = struct('Ea', Ea, 'QF', QF, 'R', R, 'QS', QS, 'QT', QT, ...
             'Sf', Sf, 'Su', Su, 'Ss', Ss, 'St', St, 'AL', AL, 'IE', IE, ...
             'SE', SE, 'Ei', Ei, 'Et', Et, 'S_canopy', S_canopy, ...
             'pot_inf', pot_inf, 'ET_vars', ET_vars, 'P', P, 'Ep', Ep);
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fast Flow Reservoir (Sf)

function [Qf,S_dt,ie]= SF_eq(S0,P,alpha,Tf,SE)
 
 ie    = (P*alpha);
 S_dt  = S0 + (P*alpha) + SE - (S0/Tf);
 Qf    = S0/Tf;
end



%% Saturated Zone Reservoir (Ss)

function [Qs,S_dt]= SS_eq(S0,R,Ts)
 
 S_dt = S0 + R - (S0/Ts);
 Qs   = S0/Ts;
end

%% UnSaturated Zone Reservoir (Ss)
function [r,se,S_dt,   S_canopy,E_int_mm,E_tran_mm,E_total_mm,rain_pass,r_a,r_s]  = SU_eq_mod(Su_0,Precip,Sumax,alpha,beta,  vars,   Sin,Lnet,Temp,u2,e_actual,S_canopy_old)


%% Define some functions
e_sat_calc      = @(Temp)                     0.6108 * exp((17.27 * Temp) ./ (237.3 + Temp));
ra_calc         = @(u, zm, d, z0, zm_p)       1/(0.4^2*u) * log((zm - d)/z0) * log((zm_p - d)/(z0/10));
g_r_calc        = @(Sin, KR)                  max(Sin * (1000 + KR) ./ (1000 * (Sin + KR)), 0);
g_D_calc        = @(KD1, KD2, VPD)            max(1 + KD1 * VPD + KD2 * VPD.^2, 0);

alpha_T         = @(TH, T0, TL)                     (TH - T0) / (T0 - TL);
g_T_calc        = @(Temp, TL, TH, T0, alpha_T)      max(  (Temp - TL) .* (TH - Temp).^alpha_T ./ ((T0 - TL) .* (TH - T0).^alpha_T),  0  );
% g_SM_calc       = @(KM1, KM2, SM, SMo)              max(1 - KM1 * exp(KM2 * (SM - SMo)), 0);
g_SM_calc      = @(FC, WP, SM)   min( max( (SM-WP)./(FC-WP), 0), 1);

Penman_Monteith_calc = @(Delta, A, rho_a, c_p, VPD, ra, gamma, rs)          ( Delta*A + rho_a*c_p*VPD/ra ) / ( Delta + gamma*(1 + rs/ra) );

Lambda_calc      = @(Temp)                                           10^3 * (2.501 - 0.002361 * Temp);          % Latent Heat of Vaporization (kJ/kg)
Delta_calc       = @(Temp, e_sat)                                    4098 * e_sat ./ (237.3 + Temp).^2;         % Slope of esat versus temp curve (kPa/C)
gamma_calc       = @(cp, Press, lambda)                              cp .* (Press)' ./ (0.622 .* lambda .* 10^3);
%% unpack some constants
zm      = vars.zm;
z0      = vars.z0;
d       = vars.d;

g_0     = vars.g0;
KR      = vars.KR;
KD1     = vars.KD1;
KD2     = vars.KD2;

TH      = vars.TH;
T0      = vars.T0;
TL      = vars.TL;
alpha_gt   = alpha_T(TH, T0, TL);

KM1     = vars.KM1;
KM2     = vars.KM2;

FC     = vars.FC;
WP     = vars.WP;

albedo   = vars.a;
cp       = vars.cp;
rho_air  = vars.ra;   

S_canopy_max = vars.Scanmax;
Patm         = vars.Patm;

%% Canopy Water balance
S_canopy_temp = S_canopy_old + Precip;

if S_canopy_temp > S_canopy_max
    rain_pass     = S_canopy_temp - S_canopy_max;
    S_canopy_temp = S_canopy_max;

else
    rain_pass     = 0;
end
% fprintf('S_canopy_old: %f\n', S_canopy_old);
% fprintf('Precip: %f\n', Precip);
% fprintf('S_canopy_temp: %f\n', S_canopy_temp);

%% Soil Water balance

 Su_0 = Su_0 + (rain_pass*(1-alpha));
 
 
 if Su_0<Sumax
      R = 0;
 else
      R = Su_0 - Sumax;
 end

 Su_0 = Su_0 - R;

%% ET
    r_a        = ra_calc(u2,zm,d,z0,zm);
    e_sat      = e_sat_calc(Temp);
    VPD        = e_sat - e_actual;
    Temp_K     = Temp + 273.17;

% Compute g's
g_r       = g_r_calc(Sin,KR);
g_D       = g_D_calc(KD1,KD2,VPD);
g_T        = g_T_calc(Temp_K, TL, TH, T0, alpha_gt);
% g_SM      = g_SM_calc(KM1,KM2,Su_0,Sumax);
g_SM      = g_SM_calc(FC,WP,Su_0);

g_s      = g_0*g_r.*g_D.*g_T.*g_SM;

% Compute rs
if  g_s ==0 
    r_s  = 10^6;
else
    r_s = 1000./g_s;
end

% Energy Balance
Rn       = Sin.*(1 - albedo) + Lnet;
Lambda   = Lambda_calc(Temp);
Delta    = Delta_calc(Temp, e_sat);
Gamma    = gamma_calc(cp,Patm,Lambda);
Wm2_mm   = (Lambda.*10^3).^-1.*(60*60*24);


% Canopy water balance


A             = Rn;
E_pot         = Penman_Monteith_calc(Delta, A, rho_air, cp, VPD, r_a, Gamma, 0);
E_tran        = Penman_Monteith_calc(Delta, A, rho_air, cp, VPD, r_a, Gamma, r_s);
E_pot_mm      = E_pot.*Wm2_mm;

if E_pot_mm > S_canopy_temp
    excess_ET     = E_pot_mm - S_canopy_temp;
    E_int_mm      = S_canopy_temp;
else
    excess_ET     = 0;
    E_int_mm      = E_pot_mm;
end
E_tran_mm     = E_tran.*Wm2_mm;



S_canopy      = S_canopy_temp -  E_int_mm; 
Su_0          = Su_0 - E_tran_mm;
E_total_mm    = E_int_mm + E_tran_mm;

% fprintf('Potential ET: %f\n', E_pot_mm);
% fprintf('Excess ET: %f\n', excess_ET);
% fprintf('Actual Interception: %f\n', E_int_mm);
% fprintf('Surface Resistance: %f\n', r_s);
% fprintf('Surface Conductance: %f\n', g_s);
% fprintf('gR: %f\n', g_r);
% fprintf('gD: %f\n', g_D);
% fprintf('gSM: %f\n', g_SM);
% fprintf('gT: %f\n', g_T); 
% fprintf('Actual Transpiration: %f\n', E_tran_mm);
% fprintf('S_canopy_temp - Inter: %f\n', S_canopy);
% fprintf('Wm2 to mm: %f\n\n', Wm2_mm);

 %%
 
 S_dt = Su_0; 
 r  = R*beta;
 se = R*(1-beta);

end

