function variables = water_energy_variables(option)


if option == 1
    % Forest
    variables.LAI = 4.00;         % (none)
    variables.h = 20.00;          % (m)
    variables.a = 0.12;           % (none)
    variables.g0 = 15.00;         % (mm s^-1)
    variables.SMo = 80.00;        % mm
    variables.SMinit = 40.00;     % mm
    variables.Scanmax = 4.00;           % (mm)
    variables.d = 14.644;         % (m)
    variables.z0 = 1.607;         % (m)
    variables.P = 101.20;         % (k Pa)
    variables.zm = 22.00;         % (m)
    variables.k = 0.4;            % (none)
    variables.esurface = 0.95;    % (none)
    variables.s = 5.67E-08;       % (W m^-2 K^-4)
    variables.KR = 200;           % (W m^-2)
    variables.KD1 = -0.307;       % (kPa^-1)
    variables.KD2 = 0.019;        % (kPa^-2)
    variables.TL = 273.00;        % (K)
    variables.T0 = 293.00;        % (K)
    variables.TH = 313.00;        % (K)
    variables.KM1 = 3.36E-04;     % (none)
    variables.KM2 = -0.10;        % (mm^-1)
    variables.ra = 1.23;          % (kg m^-3)
    variables.cp = 1013.00;       % J kg^-1 K^-1
    variables.gc = 1.00;          % (none)
    variables.aT = 1.00;          % (none)
end



if option == 2
    
    %Grass
    variables.LAI = 2.00;         % (none)
    variables.h = 0.12;           % (m)
    variables.a = 0.23;           % (none)
    variables.g0 = 30.00;         % (mm s^-1)
    variables.SMo = 40.00;        % mm
    variables.SMinit = 20.00;     % mm
    variables.Scanmax = 2.00;           % (mm)
    variables.d = 0.077;          % (m)
    variables.z0 = 0.013;         % (m)
    variables.P = 101.2;          % (k Pa)
    variables.zm = 2.00;          % (m)
    variables.k = 0.4;            % (none)
    variables.esurface = 0.95;    % (none)
    variables.s = 5.67E-08;       % (W m^-2 K^-4)
    variables.KR = 200;           % (W m^-2)
    variables.KD1 = -0.307;       % (kPa^-1)
    variables.KD2 = 0.019;        % (kPa^-2)
    variables.TL = 273.00;        % (K)
    variables.T0 = 293.00;        % (K)
    variables.TH = 313.00;        % (K)
    variables.KM1 = 1.87E-02;     % (none)
    variables.KM2 = -1.00E-01;    % (mm^-1)
    variables.ra = 1.23;          % (kg m^-3)
    variables.cp = 1013.00;       % J kg^-1 K^-1
    variables.gc = 1.00;          % (none)
    variables.aT = 1.00;          % (none)

end



end