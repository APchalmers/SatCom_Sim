%********************************************************************
% Wireless, Photonics and Space Engineering MSc. Programme
%
% Alvaro Perez Ortega
%
% Chalmers University
%********************************************************************
% This work is licensed under a Creative Commons 
% Attribution-ShareAlike 4.0 International License.
% Info: https://creativecommons.org/licenses/by-sa/4.0/
%********************************************************************

close all;
clear all;
clc;

%% Constants
To = 290;                   % Standard Temperature [K]
Tm = 275;                   % Temperature of the rain [K]
c0 = 3e8;                   % Speed of light [m/s]
f = 4e9;                    % Frequency[GHz]
seff = 2;                   % Polarization (Horizontal or Vertical)
tau = 45;                   % Polarization Angle
p = 0.001;                  % Acumulated rain probability
k = 1.38064852e-23;         % Boltzman constant[J/K]
lambda = c0/f;              % Wavelength [m]
N = 165*2;                  % Number of Clients per boat
Rb = 1*1e6*N;               % Bit Rate [bps]
Ebno = 9.6;                 % Energy per bit for QPSK Pe=1e-5
B = Rb/2;                   % Required Bandwidth [Hz]

%% Functions
dB = @(x)10*log10(x);
nat =@(x)10^(x/10);
F2T =@(F)(F - 1)*To;
T2F =@(Te) 1 + Te/To;
D_ant = @(theta3dB) 1.2*lambda./deg2rad(theta3dB);
A_EFF = @(D,eta) pi*D.^2/4*eta;
G = @(D,eta,lambda) 4*pi*A_EFF(D,eta)./lambda^2;

%% Main 

% Antennas
theta3dB = [10.5 5.5 8.5];      % 3dB Antenna Angles [Degrees]
eta = 0.7;                      % Antenna Efficiency [unitless]
Dant = D_ant(theta3dB);         % Antenna Diameters [m]
Gant = dB(G(Dant,eta,lambda));  % Antenna Gains [dB]

% Transmiter
Ptx = 330;                  % Transmitter Power [W]
Lftx_dB = 0.5;              % Transmitter Feeder Loss [dB]
Gtx_dB = Gant(2);           % Gain of Transmitter Antenna[dBi]
Lrtx_dB = 3;                % Missalignment Loss Transmitter[dB]

% Receiver
Tsky = 8;                       % Sky Temperature [K]
Tground = 28;                   % Ground Temperature[K]
Drx = 3;                        % Diameter Receiver Antenna [m]
Grx_dB = dB(G(Drx,eta,lambda)); % Gain of Receiver Antenna[dBi]
Lrrx_dB = 1;                    % Missalignment Loss Receiver[dB]
Tfrx = 280;                     % Physical Temperature of the feeder[K]
Lfrx_dB = 0.5;                  % Receiver Feeder Loss [dB]
Te_rx = 65;                     % Equivalent Temperature of the Receiver[K]

% Rain Attenuation
%   Nicaragua
LatT = 13;                  % Latitude of a receiver terminal [Degrees]
LonT = -84;                 % Lontitude of a receiver terminal [Degrees]
zone = 'P';                 % Rain Zone (ITU Maps)
% London
% LatT = 51.45;               % Latitude of a receiver terminal [Degrees]
% LonT = 0.36;                % Lontitude of a receiver terminal [Degrees]
% zone = 'F';                 % Rain Zone (ITU Maps)
LonS = -45;                 % Longitude of the subsatelite point [Degrees]                 
h_ant = 0;                  % Height of the receiving antenna [m]
h_sat = 35786e3;            % Satellite height (Geostationary) [m]

% Rain function - Includes distance and angles for the Ground Station
% [d,azimuth,elevation, gammaR, Arain_dB,K1,alpha1] = rainAttenuation(p,zone...
%     ,f*1e-9,tau,h_ant,h_sat,LatT,LonT,LonS);
% [d,azimuth,elevation, gammaR, Arain_dB,K2,alpha2] = rain_ITU(p,zone,...
%     f*1e-9,tau,h_ant,h_sat,LatT,LonT,LonS);
[d,azimuth,elevation, gammaR01, Arain_dB,K,alpha] = rain_ITU_2015(p,zone,...
    f*1e-9,tau,h_ant,h_sat,LatT,LonT,LonS);

% Free Space Loss [dB]
Lfs_dB = 20*log10((4*pi*d)/lambda);

% Noise Temp. of the rain [K]
Train = Tm*(1 - 1/nat(Arain_dB));  
% Antenna Temperature [K]
Ta = Tsky/nat(Arain_dB) + Train + Tground;
% System Temperature [K]
Tsys = Ta/nat(Lfrx_dB) + Tfrx*(1 - 1/nat(Lfrx_dB)) + Te_rx;

% EIRP
EIRP = dB(Ptx) - Lftx_dB - Lrtx_dB + Gtx_dB;
% GT
GT = Grx_dB - Lrrx_dB - Lfrx_dB - 10*log10(k) - 10*log10(Tsys);
% CNo
CNo = EIRP - Lfs_dB - Arain_dB + GT;
% SNR
SNR = EIRP - Lfs_dB - Arain_dB + GT - 10*log10(B);

fprintf('*************************************************************\n');
fprintf('    SYSTEM INFO - DOWNLINK - TX (Satellite) - RX (Earth)   \n\n');
fprintf('Frequency: %.2f GHz \n', f*1e-9);
fprintf('Polarization Angle: %d Degrees\n',tau);
fprintf('Clients: %d \n',N);
fprintf('Bit Rate: %.2f Mb/s \n',Rb*1e-6);
fprintf('Bandwidth: %.2f MHz \n',B*1e-6);
fprintf('Modulation: QPSK \n');
fprintf('Spectral efficiency: %d \n',seff);
fprintf('Eb/no required: %.2f dB \n',Ebno);
fprintf('SNR required: %.2f dB \n',dB(nat(Ebno)*seff));
fprintf('\n                         TRANSMITER                   \n\n');
fprintf('Transmitted Power: %.2f W | %.2f dB \n', Ptx, dB(Ptx));
fprintf('Feeder Loss = %.2f dB \n',Lftx_dB);
fprintf('Missalignment Loss = %.2f dB \n',Lrtx_dB);
fprintf('Antenna 3dB Angles = %.2f Degrees \n',theta3dB);
fprintf('Antenna Gains = %.2f dB \n',Gant);
fprintf('Antenna Diameters = %.2f m \n',Dant);
fprintf('\n                          RECEIVER                   \n\n');
fprintf('Feeder Loss = %.2f dB \n',Lfrx_dB);
fprintf('Missalignment Loss = %.2f dB \n',Lrrx_dB);
fprintf('Antenna Temperature = %.2f K \n',Ta);
fprintf('System Temperature = %.2f K \n',Tsys);
fprintf('Antenna Gain = %.2f dB \n',Grx_dB);
fprintf('Antenna Diameter = %.2f m \n',Drx);
fprintf('\n                         LINK BUDGET                   \n\n');
fprintf('Free Space Path Loss = %.2f dB \n',Lfs_dB);
fprintf('EIRP = %.2f dB \n',EIRP);
fprintf('G/T = %.2f dB \n',GT);
fprintf('C/No = %.2f dB \n',CNo);
fprintf('SNR = %.2f dB \n',SNR);
fprintf('*************************************************************\n');