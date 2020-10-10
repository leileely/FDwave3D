function [M] = FPgenMT(strike, dip, rake, gamma, sigma)
% This function calculates moment tensor using fault plane parameters and shear-tensile source model.

%   Input parameters:
%
%     strike, dip, rake: fault plane parameters (degrees).
%     gamma:  tensile angle in degrees (0 degrees for pure shear faulting, 
%             90 degrees for pure tensile opening).
%     sigma:  Poisson's ratio.

%   Output parameters:

%   References:
%
%     [1] Kwiatek, G. and Y. Ben-Zion (2013). Assessment of P and S wave 
%         energy radiated from very small shear-tensile seismic events in 
%         a deep South African mine. J. Geophys. Res. 118, 3630-3641, 
%         doi: 10.1002/jgrb.50274
%     [2] Ou, G.-B., 2008, Seismological Studies for Tensile Faults. 
%         Terrestrial, Atmospheric and Oceanic Sciences 19, 463.
%     [3] Vavrycuk, V., 2001. Inversion for parameters of tensile 
%         earthquakes.J. Geophys. Res. 106 (B8): 16339-6355. 
%         doi: 10.1029/2001JB000372.


strike = strike * pi/180;
dip = dip * pi / 180;
rake = rake * pi / 180;
gamma = gamma * pi / 180;
M=zeros(3,3);

M(3,3)=(sin(gamma).*(2.*cos(dip).^2 - (2.*sigma)./(2.*sigma - 1)) + cos(gamma).*sin(2.*dip).*sin(rake));
M(1,3)=(cos(gamma).*(-cos(2.*dip).*sin(rake).*sin(strike) - cos(dip).*cos(rake).*cos(strike)) + sin(gamma).*sin(2.*dip).*sin(strike));
M(2,3)=(cos(gamma).*(cos(2.*dip).*cos(strike).*sin(rake) - cos(dip).*cos(rake).*sin(strike)) - sin(gamma).*sin(2.*dip).*cos(strike));
M(1,2)=(cos(gamma).*(sin(dip).*cos(2.*strike).*cos(rake) + sin(2.*dip).*sin(rake).*sin(2.*strike)./2) - sin(gamma).*sin(dip).^2.*sin(2.*strike));
M(2,2)=(cos(gamma).*(sin(2.*strike).*cos(rake).*sin(dip) - (sin(2.*dip).*sin(rake).*cos(strike)).*2) - sin(gamma).*((2.*sigma)./(2.*sigma - 1) - 2.*cos(strike).^2.*sin(dip).^2));
M(1,1)=(cos(gamma).*(-sin(2.*strike).*cos(rake).*sin(dip) - sin(2.*dip).*cos(strike).^2.*sin(rake)) - sin(gamma).*((2.*sigma)./(2.*sigma - 1) - 2.*sin(strike).^2.*sin(dip).^2));
M(3,2)=M(2,3);
M(3,1)=M(1,3);
M(2,1)=M(1,2);

    
end
