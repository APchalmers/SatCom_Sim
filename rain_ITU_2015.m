function [d,azimuth,elevation, gammaR01, Arain,K,alpha] = rain_ITU_2015(p,...
    zone,freq,tau,h_ant,h_sat,LatT,LonT,LonS)
%RAINATTENUATION Summary of this function goes here
%   Detailed explanation goes here

% Elevation Calculation
[d,azimuth,elevation] = calculateCoord(LatT,LonT,0,LonS,h_sat);

% Frequencies and Ks and Alphas (ITU 838)
[K,alpha] = Kcalc(freq,tau,elevation);

% Zones and Rainfall rate, R01 (ITU 837)
zoneV = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'J' 'K' 'L' 'M' 'N' 'P' 'Q'];
R01v = [8 12 15 19 22 28 30 32 35 42 60 63 95 145 115];

% 1) Height of the rain (ITU 839)
if (LatT > 23 && LonT > -60)
    h_o = 5 - 0.075*(LatT - 23);
elseif (LatT > 0 && LatT < 23)
    h_o = 5;
elseif (LatT < 0 && LatT > -21)
    h_o = 5;
elseif (LatT < -21 && LatT > -71)
    h_o = 5 + 0.1*(LatT + 21);
elseif (LatT < -71)
    h_o = 0;
end
if (LatT > 35 && LatT < 70 && LonT < -60)
    h_rain = 3.2 - 0.075*(LatT - 35);
else
    h_rain = h_o + 0.36;
end

% 2) Distance of the signal through the rain ITU (2015)
if elevation < 5
    Ls = 2*(h_rain - h_ant)/((sind(elevation).^2 + 2*(h_rain - h_ant)...
        /earthRadius('meter')).^(1/2) + sind(elevation));
else
    Ls = (h_rain - h_ant)/sind(elevation);
end

% 3) Calculate Horizontal Projection
Lg = Ls*cosd(elevation);

% 4) Search R01 according to zone
[~,indZ] = find(zoneV == zone);
R01 = R01v(indZ);

% 5) Get specific attenuation gammaR 
gammaR01 = K*R01^alpha;

% 6) Calculate the horizontal reduction factor
r01 = 1/(1 + 0.78*sqrt(Lg*gammaR01/freq) - 0.38*(1-exp(-2*Lg)));
% 7) Calculate the vertical adjustment factor
psi = atand((h_rain - h_ant)/(Lg*r01));
if psi > elevation
    Lr = (Lg*r01)/cosd(elevation);
else
    Lr = (h_rain - h_ant)/sind(elevation);
end
if abs(LatT) < 36
    xsi = 36 - abs(LatT);
else 
    xsi =0;
end
v01 = 1/(1 + sqrt(sind(elevation))*(31*(1 - exp(-(elevation/(1 + xsi))))*...
    sqrt(Lr*gammaR01)/freq^2 - 0.45));

% 8) The effective path length
Le = Lr*v01;
% 9) The predicted rain attenuation exceeded for 0.01% of an average year
A01 = gammaR01*Le;
% 10) Scale

if (p >= 1 || abs(LatT)>36)
    beta = 0;
elseif (p < 1 && abs(LatT) < 36 && elevation >= 25)
    beta = -0.005*(abs(LatT)-36);
else
    beta = -0.005*(abs(LatT)-36) + 1.8 - 4.25*sind(elevation);
end

Arain = A01*(p/0.01)^(-(0.655+0.033*log(p) - 0.045*log(A01) - beta*...
    (1 - p)*sind(elevation)));


fprintf('*************************************************************\n');
fprintf('                         RAIN MODEL INFO                   \n\n');
fprintf('Point Coordinates: [%.2fN, %.2fE] \n',LatT,LonT);
fprintf('Rain Zone: %s \n', zoneV(indZ));
fprintf('Frequency for model: %.3f GHz \n', freq);
fprintf('Rain Rate (0.01%%) = %d mm/h \n',R01);
fprintf('Elevation = %.2f Degrees \n',elevation);
fprintf('Azimuth = %.2f Degrees \n',azimuth);
fprintf('Point-to-Point Distance = %.3f km \n',d*1e-3);
fprintf('Rain Height = %.3f km \n',h_rain);
fprintf('Slant Path = %.3f km \n',Ls);
fprintf('Effective Path = %.3f km \n',Le);
fprintf('Attenuation (gamma)(0.01%%) = %.5f [dB/km] \n',gammaR01);
fprintf('Total Attenuation (0.01%%) = %.3f [dB] \n',A01);
fprintf('Total Attenuation (scaled)(%.3f %%) = %.3f [dB] \n',p,Arain);
fprintf('*************************************************************\n');
end

