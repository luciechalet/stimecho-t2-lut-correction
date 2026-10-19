% PHANTOM_TESTS  Tests on the GRAPPATINI phantom (2026-10-12, Cima.X, syngo XA61).
%
%   Run from the repository root:   matlab -batch phantom_tests
%
%   Needs the phantom DICOM and BIDS folders, looked for in the sibling folder
%   ../test_data/GRAPPATINI/phantom. Skips, without failing, when they are absent.
%   No StimFit and no proprietary/ needed.

fprintf('\nphantom_tests: %s\n\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

here    = fileparts(mfilename('fullpath'));
phantom = fullfile(fileparts(here), 'test_data', 'GRAPPATINI', 'phantom');
dcm     = fullfile(phantom, 'DICOM');
anat    = fullfile(phantom, 'bids', 'sub-phantom', 'ses-20261012', 'anat');
if ~isfolder(dcm) || ~isfolder(anat)
    fprintf('phantom_tests skipped: phantom data not found under %s\n\n', phantom);
    return
end

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

function s = list_if_any(c)
    if isempty(c), s = ''; else, s = [' - differs: ' strjoin(c, ', ')]; end
end
