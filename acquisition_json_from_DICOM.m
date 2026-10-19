function acq = acquisition_json_from_DICOM(dicom_dir, args)
    arguments
        dicom_dir (1,:) char
        args.SeriesNumber double = []        % required when the folder holds several series
        args.OutputFile (1,:) char = ''      % '' = return and print only, write nothing
        args.verbose (1,1) logical = true
    end
    % Builds the acquisition-parameters JSON (README "Inputs" 1) from a Siemens DICOM
    % series, as populate_json_sidecar_from_PDF.py does from a protocol printout - but
    % from the DICOM alone, no sidecar needed.
    % Usage:
    %       acq = acquisition_json_from_DICOM(dicom_dir, 'SeriesNumber', 7);
    %       acq = acquisition_json_from_DICOM(dicom_dir, 'SeriesNumber', 7, 'OutputFile', 'acq.json');
    % Inputs:
    %       dicom_dir: folder holding the DICOM files (searched recursively, may be flat)
    %       args.SeriesNumber: series to read; may be omitted only if one series in the
    %                          folder carries a protocol block
    %       args.OutputFile: path of the JSON to write (default: write nothing)
    %       args.verbose: print what was read and from where (default: true)
    % Outputs:
    %       acq: struct ready for build_opt_from_protocol, with identity fields written
    %            as dcm2niix writes them so check_LUT_compatibility can compare a LUT built
    %            from it with ordinary dcm2niix sidecars
    % Dependencies: Image Processing Toolbox (dicominfo) for the standard DICOM fields

    %% Find the files of the series that carry the protocol block
    files = dir(fullfile(dicom_dir, '**', '*'));
    files = files(~[files.isdir]);
    if isempty(files)
        error('acquisition_json_from_DICOM:NoFiles', 'No files found in %s.', dicom_dir);
    end

    % The scan reads each file's bytes once and looks for two things: a protocol block,
    % and the series number. dicominfo is NOT used here: on multi-frame files (MP2RAGE,
    % SWI) it decodes every per-frame header and took minutes per folder. It runs once,
    % below, on the chosen file only, and confirms the series number read here.
    cand = struct('file', {}, 'series', {}, 'description', {});
    for k = 1:numel(files)
        f = fullfile(files(k).folder, files(k).name);
        bytes = read_bytes(f);
        if ~contains(bytes, '### ASCCONV BEGIN')
            continue    % not an image with a protocol (e.g. a presentation state), or not DICOM
        end
        sn = series_number_from_bytes(bytes);
        if isnan(sn)
            % Not explicit VR little endian: fall back to the slow, general reader
            info = dicominfo(f, 'UseDictionaryVR', true);
            sn = double(info.SeriesNumber);
        end
        name = regexp(bytes, 'tProtocolName\s*=\s*"*([^"\n]*)', 'tokens', 'once');
        if isempty(name), name = {''}; end
        cand(end+1) = struct('file', f, 'series', sn, 'description', name{1}); %#ok<AGROW>
    end
    if isempty(cand)
        error('acquisition_json_from_DICOM:NoProtocol', ...
            'No file in %s carries a Siemens protocol (ASCCONV) block.', dicom_dir);
    end

    series = unique([cand.series]);
    if isempty(args.SeriesNumber)
        if numel(series) > 1
            error('acquisition_json_from_DICOM:WhichSeries', ...
                ['%s holds %d series with a protocol block; choose one with ''SeriesNumber''.\n%s'], ...
                dicom_dir, numel(series), series_list(cand));
        end
        args.SeriesNumber = series;
    end
    sel = cand([cand.series] == args.SeriesNumber);
    if isempty(sel)
        error('acquisition_json_from_DICOM:SeriesNotFound', ...
            'Series %d has no file with a protocol block in %s.\n%s', ...
            args.SeriesNumber, dicom_dir, series_list(cand));
    end
    src = sel(1).file;    % every file of a series carries the same protocol

    %% Protocol facts from the ASCCONV block
    facts = siemens_protocol_facts(read_siemens_ascconv(src));

    %% Identity from the standard DICOM header, written as dcm2niix writes it
    % 'UseDictionaryVR' is faster on XA61 Enhanced MR files and avoids a warning about
    % a mismatched VR in a tag we do not use (CompletionFlag, (0040,A491)).
    info = dicominfo(src, 'UseDictionaryVR', true);
    if double(info.SeriesNumber) ~= args.SeriesNumber
        error('acquisition_json_from_DICOM:SeriesMismatch', ...
            'Byte scan placed %s in series %d, but its header says %d.', ...
            src, args.SeriesNumber, double(info.SeriesNumber));
    end
    if ~startsWith(info.Manufacturer, 'Siemens', 'IgnoreCase', true)
        error('acquisition_json_from_DICOM:NotSiemens', ...
            'Manufacturer is "%s"; the ASCCONV block is Siemens-specific.', info.Manufacturer);
    end
    acq.ProtocolName           = char(info.ProtocolName);
    acq.SeriesDescription      = char(info.SeriesDescription);
    acq.SeriesNumber           = double(info.SeriesNumber);
    acq.Manufacturer           = 'Siemens';          % dcm2niix shortens "Siemens Healthineers"
    acq.ManufacturersModelName = char(info.ManufacturerModelName);
    acq.MagneticFieldStrength  = double(info.MagneticFieldStrength);
    acq.SoftwareVersions       = char(info.SoftwareVersions);

    %% Acquisition fields from the protocol
    acq.PulseSequenceDetails = facts.PulseSequenceDetails;
    acq.EchoTime             = facts.EchoTime;
    acq.EchoTrainLength      = facts.EchoTrainLength;
    acq.RepetitionTime       = facts.RepetitionTime;
    acq.FlipAngle            = facts.FlipAngle;
    acq.SliceThickness       = facts.SliceThickness;
    acq.RFPulseTypeCode      = facts.RFPulseTypeCode;

    % Provenance, mirroring the SiemensProtocolPrintout block of the PDF tool
    acq.SiemensProtocolHeader = struct( ...
        'SourceFile',   src, ...
        'SeriesNumber', args.SeriesNumber, ...
        'FieldSources', facts.Source);

    %% Report
    if args.verbose
        fprintf('Series %d (%s), read from %s\n', acq.SeriesNumber, acq.SeriesDescription, src);
        fprintf('  EchoTime        %d echoes, %.4g to %.4g ms, spacing %.4g ms\n', numel(acq.EchoTime), ...
            acq.EchoTime(1)*1e3, acq.EchoTime(end)*1e3, (acq.EchoTime(2)-acq.EchoTime(1))*1e3);
        fprintf('  RepetitionTime  %g s\n', acq.RepetitionTime);
        fprintf('  FlipAngle       %g deg\n', acq.FlipAngle);
        fprintf('  SliceThickness  %g mm\n', acq.SliceThickness);
        fprintf('  RFPulseTypeCode %d (scanner code, not translated)\n', acq.RFPulseTypeCode);
        fprintf('  Protocol        %s, %s, %s\n', acq.ProtocolName, acq.PulseSequenceDetails, acq.SoftwareVersions);
    end

    %% Write JSON if requested
    if ~isempty(args.OutputFile)
        fid = fopen(args.OutputFile, 'w');
        if fid < 0
            error('acquisition_json_from_DICOM:CannotWrite', 'Could not open %s for writing.', args.OutputFile);
        end
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '%s', jsonencode(acq, 'PrettyPrint', true));
        if args.verbose
            fprintf('Acquisition JSON written to %s\n', args.OutputFile);
        end
    end
end

function c = read_bytes(file)
    % Whole file as a char row (one char per byte). '' if it cannot be opened.
    c = '';
    fid = fopen(file, 'r');
    if fid < 0
        return
    end
    closer = onCleanup(@() fclose(fid));
    c = char(fread(fid, inf, '*uint8')');
end

function sn = series_number_from_bytes(c)
    % Series Number (0020,0011), VR IS, in explicit VR little endian: the tag bytes
    % 20 00 11 00, then 'IS', a 2-byte length and the number as text. The first match
    % is the top-level element - checked against pydicom on all 78 files of the
    % 2026-10-12 phantom. NaN if absent (e.g. implicit VR files).
    sn = NaN;
    k = strfind(c, [char([32 0 17 0]) 'IS']);
    if isempty(k)
        return
    end
    k = k(1);
    len = double(c(k+6)) + 256*double(c(k+7));
    if k + 7 + len > numel(c)
        return
    end
    sn = str2double(strtrim(strrep(c(k+8:k+7+len), char(0), '')));
end

function s = series_list(cand)
    [u, i] = unique([cand.series]);
    s = strjoin(arrayfun(@(k) sprintf('  series %d: %s', u(k), cand(i(k)).description), ...
        1:numel(u), 'UniformOutput', false), newline);
end
