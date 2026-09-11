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
    %               .Unverified (identity fields the acquisition does not carry)

    identity = {'ProtocolName', 'SequenceName', 'SoftwareVersions', 'PulseSequenceDetails'};

    assert(isfield(LUT,'Protocol') && isstruct(LUT.Protocol), 'check_LUT_protocol:NoProtocol', ...
        'LUT carries no Protocol structure - was it built by build_T2_LUT?');

    matched = {};
    unverified = {};
    for k = 1:numel(identity)
        f = identity{k};
        if ~isfield(LUT.Protocol, f)
            error('check_LUT_protocol:LUTMissingField', ...
                'LUT.Protocol has no %s, so the LUT cannot be matched on it.', f);
        end
        if ~isfield(acq, f) || isempty(acq.(f))
            unverified{end+1} = f; %#ok<AGROW>
            continue
        end
        a = as_text(LUT.Protocol.(f), ['LUT.Protocol.' f]);
        b = as_text(acq.(f), ['acquisition ' f]);
        if ~strcmp(a, b)
            % Note the values are passed as arguments, never pasted into the
            % format string: PulseSequenceDetails contains % and \ characters,
            % which sprintf would otherwise try to interpret.
            error('check_LUT_protocol:ProtocolMismatch', ...
                ['%s differs between the LUT and the acquisition:\n' ...
                 '    LUT         : %s\n' ...
                 '    acquisition : %s\n' ...
                 'This LUT was built for a different protocol and must not be applied here.'], ...
                f, a, b);
        end
        matched{end+1} = f; %#ok<AGROW>
    end

    if ~any(strcmp(matched, 'ProtocolName'))
        error('check_LUT_protocol:NoProtocolName', ...
            'The acquisition carries no ProtocolName, so there is nothing to match the LUT against.');
    end

    if ~isempty(unverified)
        warning('check_LUT_protocol:Unverified', ...
            ['Matched on %s. Absent from the acquisition and therefore NOT checked: %s. ' ...
             'A protocol can be revised without its name changing, so this is a weaker match.'], ...
            strjoin(matched, ', '), strjoin(unverified, ', '));
    end

    report = struct('Matched', {matched}, 'Unverified', {unverified});
end

function s = as_text(v, what)
    % The identity fields are strings in both files. Accept the shapes jsondecode
    % can produce, and refuse anything else rather than comparing it loosely.
    if isstring(v) && isscalar(v), v = char(v); end
    if iscellstr(v) && isscalar(v), v = v{1}; end
    if ~ischar(v) || (~isrow(v) && ~isempty(v))
        error('check_LUT_protocol:NotText', '%s is not a character string.', what);
    end
    s = v;
end
