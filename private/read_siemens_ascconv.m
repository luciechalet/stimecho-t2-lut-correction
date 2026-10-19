function [P, found] = read_siemens_ascconv(dicom_file)
    % Reads the Siemens protocol text block (ASCCONV) stored inside a DICOM file.
    % Usage:
    %       [P, found] = read_siemens_ascconv(dicom_file);
    % Inputs:
    %       dicom_file: path to one DICOM file
    % Outputs:
    %       P: containers.Map from protocol parameter name to its value as text,
    %          e.g. P('alTE[0]') = '13300'. Empty map if the file has no block.
    %       found: true if the file carries an ASCCONV block
    % Notes:
    %       The block is plain text between '### ASCCONV BEGIN' and '### ASCCONV END'
    %       inside a private tag (Enhanced MR, XA: (0021,1019); classic DICOM: the CSA
    %       series header). Searching the raw bytes finds it in both cases, so no DICOM
    %       toolbox is needed here. Values are kept as text; units are those of the
    %       protocol (times in microseconds, lengths in millimetres).
    % Dependencies: base MATLAB only

    P = containers.Map('KeyType', 'char', 'ValueType', 'char');
    found = false;

    fid = fopen(dicom_file, 'r');
    if fid < 0
        error('read_siemens_ascconv:CannotRead', 'Could not open %s.', dicom_file);
    end
    closer = onCleanup(@() fclose(fid));
    txt = char(fread(fid, inf, '*uint8')');

    b = strfind(txt, '### ASCCONV BEGIN');
    if isempty(b)
        return
    end
    if numel(b) > 1
        warning('read_siemens_ascconv:SeveralBlocks', ...
            '%s holds %d ASCCONV blocks; reading the first.', dicom_file, numel(b));
    end
    e = strfind(txt(b(1):end), '### ASCCONV END');
    if isempty(e)
        error('read_siemens_ascconv:Truncated', ...
            '%s: ASCCONV block starts but never ends.', dicom_file);
    end
    block = txt(b(1):b(1)+e(1)-2);
    found = true;

    % One parameter per line: name, '=', value, separated by tabs and spaces.
    tok = regexp(block, '^\s*([^\s=#]+)\s*=\s*(.*?)\s*$', 'tokens', 'lineanchors', 'dotexceptnewline');
    for k = 1:numel(tok)
        P(tok{k}{1}) = tok{k}{2};
    end
end
