% PHANTOM_TESTS  Tests on the GRAPPATINI phantom (2026-10-12, Cima.X, syngo XA61).
%
%   Run from the repository root:   matlab -batch phantom_tests
%
%   Uses the de-identified test set shipped in data/phantom (three
%   DICOM files and two dcm2niix sidecars, see its README). 

fprintf('\nphantom_tests: %s\n\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

here    = fileparts(mfilename('fullpath'));
shipped = fullfile(here, 'data', 'phantom');
if isfolder(shipped)
    dcm  = fullfile(shipped, 'dcm');
    anat = fullfile(shipped, 'sidecars');
end
if ~isfolder(dcm) || ~isfolder(anat)
    fprintf('phantom_tests skipped: no phantom data in %s\n\n', shipped);
    return
end
fprintf('data: %s\n\n', fileparts(dcm));

%% ------------------------------------------------------------------ %%
%  acquisition_json_from_DICOM
%% ------------------------------------------------------------------ %%
fprintf('acquisition_json_from_DICOM\n');

check_error(@() acquisition_json_from_DICOM(dcm, 'verbose', false), ...
    'acquisition_json_from_DICOM:WhichSeries', 'a folder with several series refuses to guess');
check_error(@() acquisition_json_from_DICOM(dcm, 'SeriesNumber', 8, 'verbose', false), ...
    'acquisition_json_from_DICOM:SeriesNotFound', ...
    'series 8 (a presentation state, no protocol block) is not mistaken for an image');

% Run-01 and run-02 differ only in echo spacing - and share ProtocolName, so the
% series number is the only thing that tells them apart.
runs = struct('series', {7, 15}, 'run', {'run-01', 'run-02'}, 'esp', {0.0133, 0.0101});

% Every identity field must equal what dcm2niix wrote, or check_LUT_compatibility
% would refuse a LUT built from this JSON against ordinary sidecars.
identity = {'ProtocolName', 'Manufacturer', 'ManufacturersModelName', 'MagneticFieldStrength', ...
            'SoftwareVersions', 'PulseSequenceDetails', 'EchoTrainLength', 'RepetitionTime', ...
            'FlipAngle', 'SliceThickness', 'SeriesNumber'};

for r = runs
    acq = acquisition_json_from_DICOM(dcm, 'SeriesNumber', r.series, 'verbose', false);
    check(numel(acq.EchoTime) == 16 && all(abs(diff(acq.EchoTime) - r.esp) < 1e-9) ...
          && abs(acq.EchoTime(1) - r.esp) < 1e-9, ...
          sprintf('%s: 16 echoes at %.1f ms spacing, the first at one spacing', r.run, r.esp*1e3));
    check(acq.RFPulseTypeCode == 4, sprintf('%s: RFPulseTypeCode is 4', r.run));

    side = jsondecode(fileread(fullfile(anat, ...
        sprintf('sub-phantom_ses-20261012_acq-GRAPPATINI_%s_T2map.json', r.run))));
    differs = identity(~cellfun(@(f) isfield(side, f) && isequal(acq.(f), side.(f)), identity));
    check(isempty(differs), sprintf('%s: all %d identity fields equal the dcm2niix sidecar%s', ...
        r.run, numel(identity), list_if_any(differs)));
    check(~isfield(acq, 'SequenceName') && ~isfield(side, 'SequenceName'), ...
        sprintf('%s: no SequenceName on either side (XA61 does not export it)', r.run));
end

% The JSON round trip: what is written is what build_opt_from_protocol will read.
out = [tempname '.json'];
acq = acquisition_json_from_DICOM(dcm, 'SeriesNumber', 15, 'OutputFile', out, 'verbose', false);
back = jsondecode(fileread(out));
delete(out);
check(isequal(back.EchoTime, acq.EchoTime) && back.RFPulseTypeCode == 4 ...
      && strcmp(back.PulseSequenceDetails, acq.PulseSequenceDetails), ...
      'the written JSON reads back unchanged');

%% ------------------------------------------------------------------ %%
%  check_LUT_compatibility with a DICOM folder
%% ------------------------------------------------------------------ %%
fprintf('\ncheck_LUT_compatibility (DICOM input)\n');

% Stand-in for a run-02 LUT: the protocol identity and echo train a LUT built from
% series 15 would record. Building a real one needs StimFit and proprietary/, and
% check_LUT_compatibility only reads these fields.
a15  = acquisition_json_from_DICOM(dcm, 'SeriesNumber', 15, 'verbose', false);
keep = {'ProtocolName', 'SoftwareVersions', 'PulseSequenceDetails', 'MagneticFieldStrength', ...
        'SliceThickness', 'RFPulseTypeCode'};
LUT02 = struct('Protocol', struct(), 'EchoSpacing_s', 0.0101, 'EchoTrainLength', 16);
for k = 1:numel(keep), LUT02.Protocol.(keep{k}) = a15.(keep{k}); end

rep = check_LUT_compatibility_quietly(LUT02, dcm, 'SeriesNumber', 15);
check(all(ismember({'EchoTime', 'EchoTrainLength', 'RFPulseTypeCode'}, rep.Matched)) ...
      && isequal(rep.Unverified, {'SequenceName'}), ...
      'a run-02 LUT matches series 15, echo train included (only SequenceName unverified, XA61)');

check_error(@() check_LUT_compatibility_quietly(LUT02, dcm, 'SeriesNumber', 7), ...
    'check_LUT_compatibility:ProtocolMismatch', ...
    'the same LUT on series 7 is refused: identical protocol, 13.3 ms against 10.1 ms');

% The trap the DICOM input closes: the run-01 dcm2niix sidecar has no echo times,
% so the wrong-run LUT passes on it - with EchoTime reported unverified.
side01 = fullfile(anat, 'sub-phantom_ses-20261012_acq-GRAPPATINI_run-01_T2map.json');
rep = check_LUT_compatibility_quietly(LUT02, side01);
check(ismember('EchoTime', rep.Unverified), ...
    'from the run-01 sidecar the mismatch cannot be seen, and EchoTime is reported unverified');

check_error(@() check_LUT_compatibility_quietly(LUT02, dcm), ...
    'acquisition_json_from_DICOM:WhichSeries', 'a DICOM folder without SeriesNumber refuses to guess');

fprintf('\nphantom_tests: all checks passed\n\n');

%% ------------------------------------------------------------------ %%
%  Helpers
%% ------------------------------------------------------------------ %%
function check(cond, name)
    if cond
        fprintf('  ok    %s\n', name);
    else
        fprintf('  FAIL  %s\n', name);
        error('phantom_tests:Failed', '%s', name);
    end
end

function check_error(f, expected_id, name)
    got = '(no error raised)';
    try
        f();
    catch ME
        got = ME.identifier;
    end
    check(strcmp(got, expected_id), sprintf('%s [%s]', name, expected_id));
end

function r = check_LUT_compatibility_quietly(varargin)
    % The Unverified warning is expected here; only the report is checked.
    orig = warning('off', 'check_LUT_compatibility:Unverified');
    restore = onCleanup(@() warning(orig)); %#ok<NASGU>
    r = check_LUT_compatibility(varargin{:});
end

function s = list_if_any(c)
    if isempty(c), s = ''; else, s = [' - differs: ' strjoin(c, ', ')]; end
end
