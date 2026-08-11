function G = calib_G(RF, target_cm)
    %  Computes the slice select gradient for a given RF pulse and target slice thickness
    %  Usage:
    %       G = calib_G(RF, target_cm)
    %  Inputs:
    %       RF: RF waveform structure with fields:
    %           .RF: RF waveform (G)
    %           .tau: pulse duration (s)
    %           .angle: prescribed nutation angle (degrees)
    %       target_cm: desired slice thickness (cm)
    %  Outputs:
    %       G: slice select gradient (G/cm)
    %       Solves for G such that the simulated FWHM equals the prescribed slice thickness.
    RF.G = 1;
    z = linspace(-5, 5, 20001); %cm
    RF = calcflipRF(RF, z, 0);
    a = abs(RF.alpha);
    [amax, ipk] = max(a);
    half = amax / 2;
    il = find(a(1:ipk) <= half, 1, 'last');
    ir = ipk -1 +find(a(ipk:end) <= half, 1, 'first');
    if isempty(il) || isempty(ir)
        error('calib_G:ProfileTooWide', ...
        'Half-maximum not reached inside the simulation window - widen z.');
    end
    zl = interp1(a([il il+1]), z([il il+1]), half); % interpolate the crossings
    zr = interp1(a([ir-1 ir]), z([ir-1 ir]), half);
    G = (zr - zl) / target_cm;
end