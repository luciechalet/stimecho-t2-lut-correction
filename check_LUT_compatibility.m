function report = check_LUT_compatibility(LUT, acq, args)
    arguments
        LUT (1,1) struct
        acq                                  % struct, JSON file path, or DICOM folder
        args.SeriesNumber double = []        % DICOM folder only: which series to check
    end
    % Check the LUT compatibility with acquired data.
    % A LUT is only valid for the protocol it was built from: applying one built
    % for a different protocol returns a smooth, plausible and wrong map, with no
    % error and nothing in the output to reveal it. This is the check that catches it.
    % Usage:
    %       check_LUT_compatibility(LUT, acq);                         % errors on any mismatch
    %       report = check_LUT_compatibility(LUT, acq);                % also returns what was compared
    %       check_LUT_compatibility(LUT, 'sub-01_T2map.json');         % a sidecar or acquisition JSON
    %       check_LUT_compatibility(LUT, 'DICOM', 'SeriesNumber', 7);  % a Siemens DICOM folder
    % Inputs:
    %       LUT: LUT structure, e.g. jsondecode(fileread('LUT/T2LUT_....json'))
    %       acq: acquisition parameters - a struct, e.g. jsondecode(fileread('sub-01_T2map.json')),
    %            the path of such a JSON file, or a DICOM folder read by acquisition_json_from_DICOM
    %       args.SeriesNumber: series to read from a DICOM folder (required when it holds several)
    % Outputs:
    %       report: struct with .Matched (fields compared and equal) and
    %               .Unverified (identity fields one of the two does not carry)
    % Notes:
    %       The echo train is compared whenever the acquisition lists every echo
    %       (EchoTime with two or more values): always from a DICOM folder, never from a
    %       dcm2niix sidecar of a GRAPPATINI map, which carries no echo times. Two
    %       protocols can share every identity field and differ only in echo spacing.

    acq = resolve_acquisition(acq, args.SeriesNumber);

    % Fields build_T2_LUT always writes: absent from the LUT means the file is
    % not one of ours.
    core  = {'ProtocolName', 'SoftwareVersions', 'PulseSequenceDetails'};
    % Fields checked whenever both sides carry them. SequenceName is here rather
    % than in core because syngo XA does not export it, so a LUT built from XA
    % data cannot carry it. The others, recorded since 2026-09-11, are what the RF
    % waveform, the slice profile and the T1 assumption rest on, so a difference
    % invalidates the table even when the name matches.
    extra = {'SequenceName', 'MagneticFieldStrength', 'SliceThickness', 'RFPulseType', 'RFPulseTypeCode'};
    % RFPulseType (printout label, e.g. 'Low SAR') and RFPulseTypeCode (DICOM number,
    % e.g. 4) record the same fact in two vocabularies. Each is compared only with
    % itself - never label against code - and one of them matching verifies the pulse.
    rf_pair = {'RFPulseType', 'RFPulseTypeCode'};

    % Relative tolerance for the numeric fields; everything else is compared exactly.
    tol = struct( ...
        'MagneticFieldStrength', 0.05, ...  % sidecars write 2.89362 for a 3 T magnet
        'SliceThickness',        1e-3);     % any real difference changes the slice profile

    assert(isfield(LUT,'Protocol') && isstruct(LUT.Protocol), 'check_LUT_compatibility:NoProtocol', ...
        'LUT carries no Protocol structure - was it built by build_T2_LUT?');

    matched = {};
    unverified = {};
    for f = [core, extra]
        name = f{1};

        if ~isfield(LUT.Protocol, name) || isempty(LUT.Protocol.(name))
            if ismember(name, core)
                error('check_LUT_compatibility:LUTMissingField', ...
                    'LUT.Protocol has no %s, so the LUT cannot be matched on it.', name);
            end
            unverified{end+1} = name; %#ok<AGROW>
            continue
        end

        [b, found] = acq_value(acq, name);
        if ~found
            unverified{end+1} = name; %#ok<AGROW>
            continue
        end
        a = LUT.Protocol.(name);

        rtol = 0;
        if isfield(tol, name), rtol = tol.(name); end
        if ~same_value(a, b, rtol, name)
            % The values are passed as arguments, never pasted into the format
            % string: PulseSequenceDetails contains % and \ characters, which
            % sprintf would otherwise try to interpret.
            error('check_LUT_compatibility:ProtocolMismatch', ...
                ['%s differs between the LUT and the acquisition:\n' ...
                 '    LUT         : %s\n' ...
                 '    acquisition : %s\n' ...
                 'This LUT was built for a different protocol and must not be applied here.'], ...
                name, describe(a), describe(b));
        end
        matched{end+1} = name; %#ok<AGROW>
    end

    if any(ismember(rf_pair, matched))
        unverified = setdiff(unverified, rf_pair, 'stable');
    end

    % Echo train: the LUT simulates echoes at esp, 2*esp, ..., etl*esp.
    for name = {'EchoTrainLength', 'EchoTime'}
        switch compare_echoes(LUT, acq, name{1})
            case 'matched',    matched{end+1}    = name{1}; %#ok<AGROW>
            case 'unverified', unverified{end+1} = name{1}; %#ok<AGROW>
        end
    end

    if ~any(strcmp(matched, 'ProtocolName'))
        error('check_LUT_compatibility:NoProtocolName', ...
            'The acquisition carries no ProtocolName, so there is nothing to match the LUT against.');
    end

    if ~isempty(unverified)
        warning('check_LUT_compatibility:Unverified', ...
            ['Matched on %s. Absent from one of the two files and therefore NOT checked: %s. ' ...
             'A protocol can be revised without its name changing, so this is a weaker match.'], ...
            strjoin(matched, ', '), strjoin(unverified, ', '));
    end

    report = struct('Matched', {matched}, 'Unverified', {unverified});
end

function acq = resolve_acquisition(acq, series)
    % Turn a JSON path or a DICOM folder into the acquisition struct.
    if isstring(acq) && isscalar(acq), acq = char(acq); end
    if isstruct(acq)
        assert(isscalar(acq), 'check_LUT_compatibility:BadAcquisition', ...
            'The acquisition struct must be a single struct, not an array.');
        assert(isempty(series), 'check_LUT_compatibility:SeriesWithoutDICOM', ...
            'SeriesNumber only applies when the acquisition is a DICOM folder.');
        return
    end
    if ~ischar(acq)
        error('check_LUT_compatibility:BadAcquisition', ...
            'The acquisition must be a struct, a JSON file path or a DICOM folder.');
    end
    if isfolder(acq)
        acq = acquisition_json_from_DICOM(acq, 'SeriesNumber', series, 'verbose', false);
    elseif isfile(acq) && endsWith(acq, '.json', 'IgnoreCase', true)
        assert(isempty(series), 'check_LUT_compatibility:SeriesWithoutDICOM', ...
            'SeriesNumber only applies when the acquisition is a DICOM folder.');
        acq = jsondecode(fileread(acq));
    else
        error('check_LUT_compatibility:BadAcquisition', ...
            '"%s" is neither a JSON file nor a folder of DICOM files.', acq);
    end
end

function status = compare_echoes(LUT, acq, name)
    % 'matched' or 'unverified'; errors on a mismatch. Echo times in seconds,
    % compared to within 1 microsecond - the resolution of the Siemens protocol.
    status = 'unverified';
    tol = 1e-6;
    switch name
        case 'EchoTrainLength'
            if ~isfield(LUT, 'EchoTrainLength') || ~isfield(acq, 'EchoTrainLength'), return, end
            a = LUT.EchoTrainLength;  b = double(acq.EchoTrainLength);
            if ~isscalar(b) || a ~= b
                error('check_LUT_compatibility:ProtocolMismatch', ...
                    ['EchoTrainLength differs between the LUT and the acquisition:\n' ...
                     '    LUT         : %d\n    acquisition : %s\n' ...
                     'This LUT was built for a different protocol and must not be applied here.'], ...
                    a, mat2str(b));
            end
        case 'EchoTime'
            if ~isfield(LUT, 'EchoSpacing_s') || ~isfield(LUT, 'EchoTrainLength') ...
                    || ~isfield(acq, 'EchoTime') || numel(acq.EchoTime) < 2
                return    % a single TE (one echo image, or none) says nothing about the train
            end
            TE = double(acq.EchoTime(:));
            if max(TE) > 1, TE = TE * 1e-3; end     % ms -> s, as build_opt_from_protocol does
            expected = (1:LUT.EchoTrainLength).' * LUT.EchoSpacing_s;
            if numel(TE) ~= numel(expected) || any(abs(TE - expected) > tol)
                error('check_LUT_compatibility:ProtocolMismatch', ...
                    ['EchoTime differs between the LUT and the acquisition:\n' ...
                     '    LUT         : %d echoes, %.4g ms spacing (first at %.4g ms)\n' ...
                     '    acquisition : %d echoes, %.4g ms spacing (first at %.4g ms)\n' ...
                     'This LUT was built for a different protocol and must not be applied here.'], ...
                    numel(expected), LUT.EchoSpacing_s*1e3, expected(1)*1e3, ...
                    numel(TE), (TE(2)-TE(1))*1e3, TE(1)*1e3);
            end
    end
    status = 'matched';
end

function [v, found] = acq_value(acq, name)
    % RFPulseType is not a BIDS field. populate_json_sidecar_from_PDF.py lifts it
    % from the Siemens printout and nests it under SiemensProtocolPrintout, so
    % look there as well before giving up.
    v = [];  found = false;
    if isfield(acq, name) && ~isempty(acq.(name))
        v = acq.(name);  found = true;  return
    end
    block = 'SiemensProtocolPrintout';
    if isfield(acq, block) && isstruct(acq.(block)) && isscalar(acq.(block)) ...
            && isfield(acq.(block), name) && ~isempty(acq.(block).(name))
        v = acq.(block).(name);  found = true;
    end
end

function tf = same_value(a, b, rtol, name)
    a = normalise(a, ['LUT.Protocol.' name]);
    b = normalise(b, ['acquisition ' name]);
    if ischar(a) ~= ischar(b)
        % One side quoted a number the other wrote bare - compare as numbers if
        % both parse, and as a plain mismatch if they do not.
        a = maybe_number(a);
        b = maybe_number(b);
    end
    if isnumeric(a) && isnumeric(b)
        tf = abs(a - b) <= rtol * max(1, abs(a));
    elseif ischar(a) && ischar(b)
        tf = strcmp(a, b);
    else
        tf = false;
    end
end

function v = normalise(v, what)
    % Accept the shapes jsondecode can produce; refuse anything else rather than
    % comparing it loosely.
    if isstring(v) && isscalar(v), v = char(v); end
    if iscellstr(v) && isscalar(v), v = v{1}; end
    ok = (ischar(v) && (isrow(v) || isempty(v))) || (isnumeric(v) && isscalar(v));
    if ~ok
        error('check_LUT_compatibility:NotComparable', ...
            '%s is neither a character string nor a single number.', what);
    end
end

function v = maybe_number(v)
    if ischar(v)
        n = str2double(v);
        if isfinite(n), v = n; end
    end
end

function s = describe(v)
    if isnumeric(v), s = num2str(v); else, s = char(v); end
end
