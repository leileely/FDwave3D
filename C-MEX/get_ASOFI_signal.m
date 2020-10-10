function signals = get_ASOFI_signal(G,SOURCE_SHAPE, tshift, fc, a)

% tshift  fc      amp  

signals = zeros(1, G.nt);

PI = pi; 

for nt = 1:G.nt

    t = nt * G.dt;
    ts = 1.0 / fc;

    switch SOURCE_SHAPE
        case 1
            % Ricker derivative (like SOFI)
            tau = PI * (t - 1.5 * ts - tshift) / (ts);
            amp = -(((1.0 - 2.0 * tau * tau) * exp(-tau * tau)));
            
        case 2
            % fumue
            if ((t < tshift) || (t > (tshift + ts)))
                amp = 0.0;
            else
                amp = ((sin(2.0 * PI * (t - tshift) * fc) - 0.5 * sin(4.0 * PI * (t - tshift) * fc)));
            end
        case 3
            disp('Source wavelet from file SOURCE_FILE')
            amp = 0;
        case 4
            % sinus raised to the power of three
            if ((t < tshift) || (t > (tshift + ts)))
                amp = 0.0;
            else
                amp = (0.75 * PI / ts) * (pow(sin(PI * (t - tshift) / ts), 3.0));
            end 
        case 5
            % Ricker wavelet (like in SAVA code)
            tau = PI * (t - 1.5 * ts - tshift) / (ts);
            amp = (((1.0 - 2.0 * tau * tau) * exp(-tau * tau)));
    end
    signals(nt) = amp * a;
end
