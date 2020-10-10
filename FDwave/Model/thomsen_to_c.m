
function elastic_parameters = thomsen_to_c(thomsen_parameters)
% This function calculates elastic parameters from Thomsen's parameters

alfa_thomsen  = thomsen_parameters(1);
beta_thomsen  = thomsen_parameters(2);
eps_thomsen   = thomsen_parameters(3);
delta_thomsen = thomsen_parameters(4);
gamma_thomsen = thomsen_parameters(5);

sigma_thomsen  = alfa_thomsen^2./beta_thomsen^2.*(eps_thomsen-delta_thomsen);

%--------------------------------------------------------------------
% WA parameters defined in Vavrycuk, V., 1997. Elastodynamic and 
% elastosdtatic Green tensors ..., Geophys. J. Int., 130, 786-800. 

eps(1) = alfa_thomsen^2.*(delta_thomsen-2*eps_thomsen);
eps(2) = alfa_thomsen^2.*(eps_thomsen-delta_thomsen);
eps(3) = beta_thomsen^2.*gamma_thomsen;
alfa   = alfa_thomsen.*(1+eps_thomsen);
beta   = beta_thomsen;

%--------------------------------------------------------------------
% calculation of elastic parameters

c33 = alfa_thomsen^2;
c44 = beta_thomsen^2;

c11 = c33.*(1+2*eps_thomsen);
c66 = c44.*(1+2*gamma_thomsen);
c13 = sqrt(2*c33.*(c33-c44).*delta_thomsen+(c33-c44).^2)-c44; % exact formula

elastic_parameters = [c11 c33 c44 c66 c13];



