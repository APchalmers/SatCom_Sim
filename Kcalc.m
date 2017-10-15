function [K,alpha] = Kcalc(freq,tau,elevation)
%KTEST Summary of this function goes here
%   Detailed explanation goes here
% 
% elevation = 0;

ajkH = [-5.3398 -0.35351 -0.23789 -0.94158];
bjkH = [-0.10008 1.26970 0.86036 0.64552];
cjkH = [1.13098 0.45400 0.15354 0.16817];
mkH = -0.18961;
ckH = 0.71147;

ajkV = [-3.80595 -3.44965 -0.39902 0.50167];
bjkV = [0.56934 -0.22911 0.73042 1.07319];
cjkV = [0.81061 0.51059 0.11899 0.27195];
mkV = -0.16398;
ckV = 0.63297;

ajaH = [-0.14318 0.29591 0.32177 -5.37610 16.1721];
bjaH = [1.82442 0.77564 0.63773 -0.96230 -3.29980];
cjaH = [-0.55187 0.19822 0.13164 1.47828 3.43990];
maH = 0.67849;
caH = -1.95537;

ajaV = [-0.07771 0.567267 -0.20238 -48.2991 48.5833];
bjaV = [2.33840 0.95545 1.14520 0.791669 0.791459];
cjaV = [-0.76284 0.54039 0.26809 0.116226 0.116479];
maV = -0.053739;
caV = 0.83433;

kH = 0;
kV = 0;
for i=1:4
    aux_KH = ajkH(i).*exp(-((log10(freq)-bjkH(i))./cjkH(i)).^2);
    kH = kH + aux_KH;
    
    aux_KV = ajkV(i).*exp(-((log10(freq)-bjkV(i))./cjkV(i)).^2);
    kV = kV + aux_KV;
end
kH = 10.^(kH  + mkH*log10(freq) + ckH);
kV = 10.^(kV  + mkV*log10(freq) + ckV);

aH = 0;
aV = 0;
for i=1:5
    aux_aH = ajaH(i).*exp(-((log10(freq)-bjaH(i))./cjaH(i)).^2);
    aH = aH + aux_aH;
    
    aux_aV = ajaV(i).*exp(-((log10(freq)-bjaV(i))./cjaV(i)).^2);
    aV = aV + aux_aV;
end
aH = aH + maH*log10(freq) + caH;
aV = aV + maV*log10(freq) + caV;

K = (kH + kV + (kH - kV).*cosd(elevation).^2 .* cosd(2*tau))./2;
alpha = (kH.*aH + kV.*aV + (kH.*aH - kV.*aV).*cosd(elevation).^2.*...
    cosd(2*tau))./(2*K); 

% figure();
% loglog(freq,kH);
% grid on;
% figure();
% loglog(freq,kV);
% grid on;
% figure();
% semilogx(freq,aH);
% set(gca,'Ylim',[0.4,1.8]);
% grid on;
% figure();
% semilogx(freq,aV);
% set(gca,'Ylim',[0.4,1.8]);
% grid on;
end

