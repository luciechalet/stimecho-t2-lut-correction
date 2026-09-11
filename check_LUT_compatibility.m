function report = check_LUT_compatibility(LUT, acq)
    arguments
        LUT (1,1) struct
        acq (1,1) struct
    end
    % Check the LUT compatibility with acquired data.
    % A LUT is only valid for the protocol it was built from: applying one built
    % for a different protocol returns a smooth, plausible and wrong map, with no
    % error and nothing in the output to reveal it. This is the check that catches it.
    % Usage:
    %       check_LUT_compatibility(LUT, acq);            % errors on any mismatch
    %       report = check_LUT_compatibility(LUT, acq);   % also returns what was compared
    % Inputs:
    %       LUT: LUT structure, e.g. jsondecode(fileread('LUT/T2LUT_....json'))
    %       acq: acquisition parameters, e.g. jsondecode(fileread('sub-01_T2map.json'))
    % Outputs:
    %       report: struct with .Matched (fields compared and equal) and
    %               .Unverified (identity fields one of the two does not carry)

    % Fields build_T2_LUT always writes: absent from the LUT means the file is
    % not one of ours.
    core  = {'ProtocolName', 'SequenceName', 'SoftwareVersions', 'PulseSequenceDetails'};
    % Fields recorded since 2026-09-11, checked whenever both sides carry them.
    % These are what the RF waveform, the slice profile and the T1 assumption
    % rest on, so a difference invalidates the table even when the name matches.
    extra = {'MagneticFieldStrength', 'SliceThickness', 'RFPulseType'};

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
