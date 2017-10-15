function [d, azimuth, elevation] = calculateElevation(LatT,LonT,LatS,...
    LonS,h_sat)
%CALCULATEELEVATION Summary of this function goes here
%   Detailed explanation goes here

% Due to the spherical laws of cosines for NP-S-T:
%   Latitude from 90 (North) to -90 (South)
%   Longitude from -90 (West) to 90 (East)
DiffLon = abs(LonS - LonT);
delta = acosd(cosd(90-LatS)*cosd(90-LatT) + sind(90-LatS)*sind(90-LatT)...
    *cosd(DiffLon));

% Again, apply spherical law of cosine to NP-S-T
%   Azimuth (From North = 0 degrees to East)
beta = acosd((cosd(90-LatS)-cosd(90-LatT)*cosd(delta))/(sind(90-LatT)*...
    sind(delta)));
if LonT > LonS
    azimuth = 360 - beta;
else 
    azimuth = beta;
end

% Apply planar law of cosine for the plane of the angle delta
%   Elevation (From 0 = parallel to earth surface to 90 = normal to earth)
%   d = Slant range = distance from terminal to satellite
Re = earthRadius('meter');
d = sqrt(Re^2 + (Re + h_sat)^2 - 2*Re*(Re + h_sat)*cosd(delta));
elevation = 90 - asind((Re + h_sat)*sind(delta)/d);

end

