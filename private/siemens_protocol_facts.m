function facts = siemens_protocol_facts(P)
    % Converts ASCCONV protocol parameters into acquisition-JSON fields, in BIDS units.
    % Usage:
    %       facts = siemens_protocol_facts(read_siemens_ascconv(dicom_file));
    % Inputs:
    %       P: containers.Map returned by read_siemens_ascconv
    % Outputs:
    %       facts: struct with EchoTime (s, column), EchoTrainLength, RepetitionTime (s),
    %              FlipAngle (deg), SliceThickness (mm), PulseSequenceDetails,
    %              RFPulseTypeCode, and Source (which parameter each field came from)
    % Dependencies: base MATLAB only

    if P.Count == 0
        error('siemens_protocol_facts:NoProtocol', 'Protocol map is empty - no ASCCONV block was read.');
    end

    %% Echo train
    % alTE is a fixed-size array: only the first lContrasts entries are echoes of this
    % acquisition. The rest are leftover template values (1700..3200 us on the XA61
    % phantom) that the scanner ignores.
    nEcho = get_num(P, 'lContrasts');
    TE = zeros(nEcho, 1);
    for k = 1:nEcho
        TE(k) = get_num(P, sprintf('alTE[%d]', k-1)) * 1e-6;   % us -> s
    end
    if any(diff(TE) <= 0)
        error('siemens_protocol_facts:EchoesNotIncreasing', ...
            'The first %d alTE values are not strictly increasing: %s', nEcho, mat2str(TE*1e3));
    end

    facts.EchoTime        = round(TE, 9);   % keep float noise out, as the PDF tool does
    facts.EchoTrainLength = nEcho;
    facts.RepetitionTime  = round(get_num(P, 'alTR[0]') * 1e-6, 9);
    facts.FlipAngle       = get_num(P, 'adFlipAngleDegree[0]');
    facts.SliceThickness  = get_num(P, 'sSliceArray.asSlice[0].dThickness');
    facts.PulseSequenceDetails = strip(get_text(P, 'tSequenceFileName'), '"');

    % Kept as the scanner's code, never translated: the DICOM does not name it, and the
    % label shown in printouts ('Low SAR') is not confirmed to map onto this number.
    facts.RFPulseTypeCode = get_num(P, 'sTXSPEC.ucRFPulseType');

    facts.Source = struct( ...
        'EchoTime',             sprintf('alTE[0..%d] (us), lContrasts = %d', nEcho-1, nEcho), ...
        'EchoTrainLength',      'lContrasts', ...
        'RepetitionTime',       'alTR[0] (us)', ...
        'FlipAngle',            'adFlipAngleDegree[0]', ...
        'SliceThickness',       'sSliceArray.asSlice[0].dThickness (mm)', ...
        'PulseSequenceDetails', 'tSequenceFileName', ...
        'RFPulseTypeCode',      'sTXSPEC.ucRFPulseType');
end

function v = get_text(P, key)
    if ~isKey(P, key)
        error('siemens_protocol_facts:MissingParameter', 'Protocol has no %s.', key);
    end
    v = P(key);
end

function v = get_num(P, key)
    s = get_text(P, key);
    if startsWith(s, '0x')
        v = hex2dec(s(3:end));   % ASCCONV writes some flags in hexadecimal
    else
        v = str2double(s);
    end
    if isnan(v)
        error('siemens_protocol_facts:NotANumber', '%s = "%s" is not a number.', key, s);
    end
end
