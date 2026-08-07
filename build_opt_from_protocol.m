function opt=build_opt_from_protocol(acq_json_path, pulse_json_path, RF_exc_path, RF_ref_path, args)
    arguments
        acq_json_path   (1,:) char
        pulse_json_path (1,:) char
        RF_exc_path     (1,:) char
        RF_ref_path     (1,:) char
        args.T1               = 1.33      % seconds; scalar or function handle of T2 % GM at 3T, Hagiwara et al. 2021
        args.verbose    (1,1) logical = false
    end
    
    % Dependencies: StimFit toolbox
    require_stimfit(); % add StimFit to path if not already there
    % Reads acquisition and pulse sequence parameters from the JSON
    % Reads RF waveforms from proprietary files
    % Returns opt structure created with StimFit optset and populated by parameters.

    %  Inputs:
    %   acq_json_path: path to JSON sidecar of the acquisition (see JSON minimum requirements in README table)
    %   pulse_json_path: path to JSON file containing pulse sequence descriptors (pulse angles and durations)
    %   RF_exc_path: path to RF excitation waveform file (provided by the manufacturer)
    %   RF_ref_path: path to RF reference waveform file (provided by the manufacturer)
    %   args.T1: T1 value in seconds (scalar) or function handle of T2 (default: 1.33 s for GM at 3T, Hagiwara et al. 2021)
    %   args.verbose: plots slice profiles and prints information if true (default: false)

    %  Outputs:
    %   opt: StimFit options structure populated with the acquisition, pulse sequence parameters and corresponding slice profiles
    
    %% Read JSON files
    AcqParams   = jsondecode(fileread(acq_json_path));
    PulseParams = jsondecode(fileread(pulse_json_path));

    %% Acquisition parameters: check for minimum requirements
    required_fields = {'EchoTime', 'FlipAngle', 'SliceThickness'};
    for req = 1:numel(required_fields)
        if ~isfield(AcqParams, required_fields{req})
            error('build_opt_from_protocol:MissingField', ...
                'Acquisition JSON must contain %s.', required_fields{req});
        end
    end

    TE = double(AcqParams.EchoTime(:)); 
    if numel(TE) < 2
        error('build_opt_from_protocol:SingleEcho', ...
            ['EchoTime must list every echo in the train, got %d value(s). ' ...
            'A sidecar written for a single echo is not enough.'], numel(TE)); 
    end

    if max(TE) > 1.0
        TE = TE * 1e-3; % convert to seconds if in ms
    end

    esp=TE(2)-TE(1); % echo spacing
    if max(abs(diff(TE)-esp)) > 1e-9
        error('build_opt_from_protocol:NonUniformEchoes', ...
            ['Echo times are not uniformly spaced (esp = %g s, largest deviation ' ...
            '%g s). FSEsig only produces echoes at esp:esp:esp*etl, so the ' ...
            'simulated train would not correspond to the acquired one.'], ...
            esp, max(abs(diff(TE) - esp)));
    end

    if abs(TE(1) - esp) > 1e-9
        error('build_opt_from_protocol:FirstEchoOffset', ...
            ['First echo is at %g s but the echo spacing is %g s. FSEsig always ' ...
            'emits its first echo at esp, so the simulated train would be offset ' ...
            'from the acquired one.'], TE(1), esp);
    end

    if isfield(AcqParams, 'EchoTrainLength') && AcqParams.EchoTrainLength ~= numel(TE)
        error('build_opt_from_protocol:EchoCountMismatch', ...
            'EchoTrainLength is %d but EchoTime lists %d echoes.', ...
            AcqParams.EchoTrainLength, numel(TE));
    end

    %% Pulse sequence parameters: check for minimum requirements
    required_RF_fields = {'Excitation.FlipAngle', 'Refocusing.FlipAngle', 'Excitation.PulseDuration', 'Refocusing.PulseDuration'};
    for req_RF = 1:numel(required_RF_fields)
        field_parts = strsplit(required_RF_fields{req_RF}, '.');
        if ~isfield(PulseParams, field_parts{1}) || ~isfield(PulseParams.(field_parts{1}), field_parts{2})
            error('build_opt_from_protocol:MissingField', ...
                'Pulse JSON must contain %s.', required_RF_fields{req_RF});
        end
    end

    tau_exc = duration_to_s(PulseParams, 'Excitation');
    tau_ref = duration_to_s(PulseParams, 'Refocusing');

    if abs(PulseParams.Refocusing.FlipAngle - AcqParams.FlipAngle) > 1e-9
        error('build_opt_from_protocol:FlipAngleMismatch', ...
            ['Refocusing flip angle disagrees: %g deg in the pulse JSON but %g deg ' ...
            'in the acquisition JSON. The wrong pulse file has been paired with ' ...
            'this protocol.'],PulseParams.Refocusing.FlipAngle, AcqParams.FlipAngle);
    end


    %% Load defaults StimFit options
    opt = StimFit_optset('s'); % 's' = slice-selective refocusing
    opt.esp   = esp;
    opt.etl   = numel(TE);
    opt.debug = 0;                     % FSEsig opens a figure on every call otherwise            
    if isa(args.T1,'function_handle')
        opt.T1 = args.T1;
    elseif isnumeric(args.T1) && isscalar(args.T1) && args.T1 > 0
        if args.T1 > 20
            error('build_opt_from_protocol:T1Units', ...
                'T1 = %g looks like milliseconds; it must be in seconds.', args.T1);
        end
        T1val  = args.T1;          % capture the value, not the struct
        opt.T1 = @(T2) T1val;
    else
        error('build_opt_from_protocol:BadT1', ...
            'T1 must be a positive scalar (s) or a function handle of T2.');
    end

    %% Define the slice profile simulation window
    d_slice=AcqParams.SliceThickness*0.1; % slice thickness in cm
    window=2*d_slice; % simulate  +/- 2 slice thickness
    step=0.016; % spatial sampling interval, may need adjustments for resolutions finer than 2mmISO
    opt.Dz=window*[-1 1]; % spatial window
    opt.Nz=2*ceil(window/step)+1; 

    %% Load RF waveforms
    RF_exc = read_waveform(RF_exc_path, opt.Nrf);
    RF_ref = read_waveform(RF_ref_path, opt.Nrf);

    %% Running bloch simulations to compute slice select Gradients
    opt.RFe.RF = scaleRF(RF_exc, tau_exc, PulseParams.Excitation.FlipAngle);
    opt.RFe.tau = tau_exc;
    opt.RFe.angle = PulseParams.Excitation.FlipAngle;
    opt.RFe.path = RF_exc_path;
    opt.RFe.G = calib_G(opt.RFe, d_slice); % no background gradient for excitation
    
    opt.RFr.RF = scaleRF(RF_ref, tau_ref, PulseParams.Refocusing.FlipAngle);
    opt.RFr.tau = tau_ref;
    opt.RFr.angle = PulseParams.Refocusing.FlipAngle;
    opt.RFr.path = RF_ref_path;
    opt.RFr.G = calib_G(opt.RFr, d_slice); % no background gradient for refocusing

    % Populate RF alpha with bloch simulation (FA at each slice position: slice profile expressed as angle)
    opt = getRF(opt);

    % Verbose output
    if args.verbose
        plot_slice_profiles(opt, AcqParams.SliceThickness);
    end

end

function tau = duration_to_s(RFp, which)
    if ~isfield(RFp, which)
        error('build_opt_from_protocol:MissingField', ...
            'Pulse JSON must contain a %s section.', which);
    end
    P = RFp.(which);
    if ~isfield(P,'FlipAngle') || ~isfield(P,'PulseDuration')
        error('build_opt_from_protocol:MissingField', ...
            '%s must have both FlipAngle and PulseDuration.', which);
    end
    if ~isfield(P,'PulseDurationUnits')
        error('build_opt_from_protocol:MissingUnits', ...
            '%s.PulseDurationUnits is required - refusing to guess the unit.', which);
    end
    switch lower(P.PulseDurationUnits)
        case {'s','sec','seconds'},        tau = P.PulseDuration;
        case {'ms','milliseconds'},        tau = P.PulseDuration * 1e-3;
        case {'us','microseconds'},        tau = P.PulseDuration * 1e-6;
        otherwise
            error('build_opt_from_protocol:BadUnits', ...
                'Unrecognised %s.PulseDurationUnits: "%s".', which, P.PulseDurationUnits);
    end
end

function RF = read_waveform(RF_path, Nrf)
    if exist(RF_path,"file") ~= 2
        error('build_opt_from_protocol:NoWaveform','RF waveform not found: %s', RF_path);
    end
    [~,~,extension] = fileparts(RF_path);
    switch lower(extension)
        case ".pta", [~,RF] = readRF_Siemens(RF_path, Nrf);
        case ".rho", [~,RF] = readRF_GE(RF_path, Nrf);
        otherwise, [~,RF] = readRF(RF_path, Nrf);
    end
end

function plot_slice_profiles(opt, slice_thickness_mm)
%   QC plot: simulated flip-angle profiles against the prescribed slice.
z  = linspace(opt.Dz(1), opt.Dz(2), opt.Nz) * 10;   % cm -> mm
ae = abs(opt.RFe.alpha) * 180/pi;                   % alpha is in radians
ar = abs(opt.RFr.alpha) * 180/pi;

figure('Name','build_opt_from_protocol: slice profiles');
plot(z, ae, 'LineWidth', 1.5); hold on;
plot(z, ar, 'LineWidth', 1.5);
yline(max(ae)/2, '--'); yline(max(ar)/2, '--');
xline(-slice_thickness_mm/2, ':'); xline(slice_thickness_mm/2, ':');
hold off; grid on;
xlabel('Position (mm)'); ylabel('Flip angle (deg)');
legend('excitation','refocusing','Location','best');
title(sprintf('Slice profiles (prescribed %.2f mm)', slice_thickness_mm));

fprintf('exc: peak %6.2f deg (nominal %3g), G = %.4f G/cm\n', ...
        max(ae), opt.RFe.angle, opt.RFe.G);
fprintf('ref: peak %6.2f deg (nominal %3g), G = %.4f G/cm\n', ...
        max(ar), opt.RFr.angle, opt.RFr.G);
end
