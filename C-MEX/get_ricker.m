function x = get_ricker(t,tau,fc)

%% Authors are temporarily removed for review.

% Ricker wavelet in time
% fc is peak frequency in Hz
% t is time array

    tcut = 1.5/fc;    
    s = (t-tau-tcut).*fc;
    x = (1-2*pi^2*s.^2).*exp(-pi^2*s.^2);
    %x(abs(s)>4.) = 0;    
end

