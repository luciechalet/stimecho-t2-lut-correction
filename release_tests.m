% RELEASE_TESTS  Public test tier for stimecho-t2-lut-correction (decision D12).
%
%   Run from the repository root:   matlab -batch release_tests
%
%   Covers everything that runs without proprietary data: apply_T2_LUT,
%   fit_monoexp, parse_skipEcho, and the LUT files shipped under LUT/.
%   build_T2_LUT / build_opt_from_protocol need StimFit and the vendor RF
%   waveforms, so they cannot be tested here at all - their regression baseline
%   is checkpoint_step4.m, which needs proprietary/.
%
%   Fails loudly on the first broken check and returns a nonzero exit code, so
%   it can be used as a CI step once the repository is public.

fprintf('\nrelease_tests: %s\n\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

%% ------------------------------------------------------------------ %%
%  parse_skipEcho - base MATLAB only
%% ------------------------------------------------------------------ %%
fprintf('parse_skipEcho\n');

check(isequal(parse_skipEcho(''),    []), 'empty string skips nothing');
check(isequal(parse_skipEcho([]),    []), 'empty array skips nothing');
check(isequal(parse_skipEcho('1'),    1), 'single index as text');
check(isequal(parse_skipEcho(1),      1), 'single index as number');
check(isequal(parse_skipEcho(' 1 , 3 '), [1 3]), 'list, with whitespace');
check(isequal(parse_skipEcho([1 3]), [1 3]), 'list as numeric vector');

% Regression: str2double('1,3') returns 13, not NaN, because MATLAB reads the
% comma as a thousands separator. The version this was lifted from tested
% str2double before splitting, so '1,3' silently skipped echo 13 - or nothing
% at all when the train was shorter. Keep this check.
check(~isequal(parse_skipEcho('1,3'), 13), 'comma is a separator, not a thousands mark');

check_error(@() parse_skipEcho('0'),   'parse_skipEcho:BadSpec', 'echo 0 rejected (indices are 1-based)');
check_error(@() parse_skipEcho('-1'),  'parse_skipEcho:BadSpec', 'negative index rejected');
check_error(@() parse_skipEcho('2.5'), 'parse_skipEcho:BadSpec', 'non-integer rejected');
check_error(@() parse_skipEcho('two'), 'parse_skipEcho:BadSpec', 'non-numeric text rejected');

%% ------------------------------------------------------------------ %%
%  fit_monoexp - needs the Optimization Toolbox (lsqnonlin)
%% ------------------------------------------------------------------ %%
fprintf('\nfit_monoexp\n');

TE_ms = (11.1:11.1:11.1*16).';        % the GRAPPATINI train, in ms
T2_true = 80;  S0_true = 120;
S_pure  = S0_true * exp(-TE_ms/T2_true);

% Both fitters minimise different objectives and agree only on a perfectly
% monoexponential curve - which makes this one input a test of both at once.
[t2_nls, s0_nls] = fit_monoexp(S_pure, TE_ms, 'method', 'nls');
[t2_log, s0_log] = fit_monoexp(S_pure, TE_ms, 'method', 'loglinear');
check(abs(t2_nls - T2_true) < 1e-6,  sprintf('nls recovers a pure exponential (%.10f ms)', t2_nls));
check(abs(t2_log - T2_true) < 1e-9,  sprintf('loglinear recovers it in closed form (%.12f ms)', t2_log));
check(abs(s0_nls - S0_true) < 1e-5 && abs(s0_log - S0_true) < 1e-8, 'both recover S0');

% Unit-agnostic: T2 comes back in whatever units TE was given in. This holds
% only because TolX is set tight enough that the relative TolFun governs -
% lsqnonlin's TolX is an ABSOLUTE step tolerance, so StimFit's TolX = 5e-4,
% chosen for T2 in seconds, drifts here (measured: 3e-3 ms on this curve).
% Two properties of the input are needed for that to show up at all, and both
% are deliberate: a bi-exponential, because on a perfect exponential the log-
% linear starting guess is already the answer and the solver never iterates;
% and FSEsig's amplitude of ~0.2, because with an amplitude of ~100 the step
% norm is dominated by A and no T2 tolerance can bite.
S_bi   = 0.2*(0.8*exp(-TE_ms/70) + 0.2*exp(-TE_ms/40));
t2_bi  = fit_monoexp(S_bi, TE_ms,      'method', 'nls');
t2_bis = fit_monoexp(S_bi, TE_ms/1000, 'method', 'nls');
check(abs(t2_bis*1000 - t2_bi) < 1e-5, sprintf('same answer fitted in s and in ms (%.2e ms apart)', abs(t2_bis*1000 - t2_bi)));

% Scale-invariant: every starting value is derived from the data.
t2_scaled = fit_monoexp(S_pure*1e4, TE_ms, 'method', 'nls');
check(abs(t2_scaled - t2_nls) < 1e-6, 'amplitude scale does not change T2');

% skipEcho must actually change which echoes are fitted. Corrupt echo 1 only:
% skipping it must recover the truth, keeping it must not.
S_bad    = S_pure;  S_bad(1) = S_bad(1)*1.30;
t2_skip  = fit_monoexp(S_bad, TE_ms, 'skipEcho', 1);
t2_noskip= fit_monoexp(S_bad, TE_ms, 'skipEcho', []);
check(abs(t2_skip - T2_true) < 1e-6, 'skipEcho 1 excludes the corrupted echo');
check(abs(t2_noskip - T2_true) > 1,  'and without it the corruption is visible');

check_error(@() fit_monoexp(S_pure, TE_ms(1:end-1)), 'fit_monoexp:SizeMismatch',   'S and TE length mismatch');
check_error(@() fit_monoexp(S_pure, TE_ms, 'skipEcho', 17), 'fit_monoexp:SkipOutOfRange', 'skipping past the end of the train');
check_error(@() fit_monoexp(S_pure(1:4), TE_ms(1:4), 'skipEcho', [1 3]), 'fit_monoexp:TooFewEchoes', 'fewer than 3 echoes left');

%% ------------------------------------------------------------------ %%
%  apply_T2_LUT - base MATLAB only
%% ------------------------------------------------------------------ %%
fprintf('\napply_T2_LUT\n');

% Synthetic table: a smooth monotone stand-in, so these checks test the code
% and not the contents of any particular LUT. The shipped LUT gets its own
% section below.
LUT = synthetic_LUT();
appAxis = LUT.ApparentT2(:);
trueAxis= LUT.TrueT2(:);

% Round trip: feeding the table its own apparent axis must return its own true
% axis. This is the check that the inversion is the right way round - swapping
% the two interp1 arguments runs cleanly and moves T2 the wrong way.
back = apply_T2_LUT(appAxis, LUT);
check(max(abs(back - trueAxis)) < 1e-9, sprintf('round trip through the table (max %.2e ms)', max(abs(back - trueAxis))));
check(all(back < appAxis), 'the correction moves T2 DOWN (apparent > true)');

% Out of range is NaN at both ends: no clipping, no extrapolation (D8).
edge = [appAxis(1) - 1e-6; appAxis(end) + 1e-6; 0.1; 1e4];
check(all(isnan(apply_T2_LUT(edge, LUT))), 'out-of-range voxels become NaN, not clipped');

% Shape in, same shape out - the function makes no assumption about volumes.
vol = reshape(linspace(appAxis(1), appAxis(end), 24), [2 3 4]);
out = apply_T2_LUT(vol, LUT);
check(isequal(size(out), size(vol)), 'any input shape is preserved');
check(~any(isnan(out(:))), 'and everything inside the table is corrected');

% NaN in, NaN out, and it is counted as neither below nor above range.
[~, iNaN] = apply_T2_LUT([80; NaN; 100], LUT);
check(iNaN.ConsideredVoxels == 2 && iNaN.BelowRangeVoxels == 0 && iNaN.AboveRangeVoxels == 0, 'NaN input is ignored, not counted out of range');

% Mask: outside is NaN, so one value means one thing - "not corrected".
m = [true false true];
[mv, iM] = apply_T2_LUT([80 80 80], LUT, 'mask', m);
check(isnan(mv(2)) && all(~isnan(mv([1 3]))), 'voxels outside the mask are NaN');
check(iM.ConsideredVoxels == 2 && iM.CorrectedVoxels == 2, 'info counts only masked-in voxels');
check_error(@() apply_T2_LUT([80 80 80], LUT, 'mask', [true false]), 'apply_T2_LUT:MaskSizeMismatch', 'mask must be on the T2 grid');

% info adds up on a deliberately mixed input.
mixed = [appAxis(1)-5; 60; 90; appAxis(end)+5];
[~, i2] = apply_T2_LUT(mixed, LUT);
check(i2.ConsideredVoxels == 4 && i2.BelowRangeVoxels == 1 && i2.AboveRangeVoxels == 1 && i2.CorrectedVoxels == 2, 'below + above + corrected accounts for every voxel');

% LUT validation.
badLUT = rmfield(LUT, 'TrueT2');
check_error(@() apply_T2_LUT(80, badLUT), 'apply_T2_LUT:MissingField', 'missing LUT field rejected');
sLUT = LUT;  sLUT.T2Units = 's';
check_error(@() apply_T2_LUT(80, sLUT), 'apply_T2_LUT:WrongUnits', 'a LUT in seconds is rejected');
fLUT = LUT;  fLUT.ApparentT2 = flip(fLUT.ApparentT2);
check_error(@() apply_T2_LUT(80, fLUT), 'apply_T2_LUT:NotInvertible', 'a non-monotone table is rejected');
lLUT = LUT;  lLUT.TrueT2 = lLUT.TrueT2(1:end-1);
check_error(@() apply_T2_LUT(80, lLUT), 'apply_T2_LUT:LengthMismatch', 'axes of different length rejected');

% The documented route from a saved LUT to a struct: jsondecode(fileread(...)).
% Decoding turns row vectors into columns and that must not matter.
tmp = [tempname '.json'];
fid = fopen(tmp, 'w');  fprintf(fid, '%s', jsonencode(LUT, 'PrettyPrint', true));  fclose(fid);
decoded = jsondecode(fileread(tmp));
delete(tmp);
check(isequal(apply_T2_LUT(appAxis, decoded), back), 'a LUT round-tripped through JSON behaves identically');

%% ------------------------------------------------------------------ %%
%  LUT/ - the tables actually shipped with the repository
%% ------------------------------------------------------------------ %%
fprintf('\nshipped LUT files\n');

here  = fileparts(mfilename('fullpath'));
files = dir(fullfile(here, 'LUT', 'T2LUT_*.json'));
check(~isempty(files), sprintf('LUT/ contains at least one table (%d found)', numel(files)));

for f = files.'
    L = jsondecode(fileread(fullfile(f.folder, f.name)));
    tag = f.name;

    for fld = {'TrueT2','ApparentT2','T2Units','Protocol','FitMethod','SkipEcho','T1_s','B1'}
        check(isfield(L, fld{1}), sprintf('%s: has %s', tag, fld{1}));
    end
    for fld = {'ProtocolName','SequenceName','SoftwareVersions','PulseSequenceDetails'}
        check(isfield(L.Protocol, fld{1}), sprintf('%s: Protocol.%s, the fields the pipeline matches on', tag, fld{1}));
    end

    t = L.TrueT2(:);  a = L.ApparentT2(:);
    check(strcmpi(L.T2Units,'ms'), sprintf('%s: T2Units is ms', tag));
    check(numel(t) == numel(a) && numel(t) > 1, sprintf('%s: axes are the same length (%d points)', tag, numel(t)));
    check(issorted(a,'strictascend'), sprintf('%s: apparent axis strictly increasing, so it inverts', tag));
    check(all(a > t), sprintf('%s: apparent > true everywhere, the sign physics predicts', tag));
    check(all(isfinite(a)) && all(isfinite(t)), sprintf('%s: no NaN or Inf in either axis', tag));

    % The table must survive its own apply step.
    check(max(abs(apply_T2_LUT(a, L) - t)) < 1e-9, sprintf('%s: round trips through apply_T2_LUT', tag));

    % Load-bearing and easy to break by widening T2grid: the GRAPPATINI export
    % saturates at the int16 ceiling of 4095 = 409.5 ms, and those voxels are
    % storage limits rather than measurements. They fall outside the table only
    % because its apparent maximum stays below 409.5 - extend the grid past
    % that and every saturated voxel silently gets a confident corrected value.
    check(a(end) < 409.5, sprintf('%s: apparent max %.2f ms stays below the 409.5 ms int16 ceiling', tag, a(end)));
end

% Probe values for the table the repository ships. Hardcoded here rather than
% recomputed: recomputing would need StimFit and proprietary/, and reading them
% back out of the same file would only compare the file to itself.
shipped = 'T2LUT_T2MAP_adapted_from_GRAPPATINI_T1-1330ms_B1-100_skip1_nls.json';
check(isfile(fullfile(here,'LUT',shipped)), 'the GRAPPATINI table is the one that ships');
G = jsondecode(fileread(fullfile(here,'LUT',shipped)));

check(strcmp(G.FitMethod,'nls') && isequal(G.SkipEcho,1) && G.B1 == 1 && abs(G.T1_s - 1.33) < 1e-12, ...
      'built with nls, echo 1 skipped, B1 = 1, T1 = 1.33 s (D1, D3, D5)');
check(isequal(G.TrueT2(:).', 10:1:300), 'grid is 10:1:300 ms, 291 points (D2)');

probes = [ 30  32.8350
           80  85.8161
          130 138.0322
          200 210.7088
          300 313.7600];
for k = 1:size(probes,1)
    got = G.ApparentT2(G.TrueT2 == probes(k,1));
    check(abs(got - probes(k,2)) < 5e-4, sprintf('true %3g ms -> apparent %.4f ms', probes(k,1), got));
end

slope = diff(G.ApparentT2(:)) ./ diff(G.TrueT2(:));
check(abs(min(slope) - 1.0261) < 5e-4, sprintf('min local slope %.4f - above 1, so the inversion never amplifies noise', min(slope)));

%% ------------------------------------------------------------------ %%
%  check_LUT_compatibility - base MATLAB only
%% ------------------------------------------------------------------ %%
fprintf('\ncheck_LUT_compatibility\n');

% The acquisition side is written out here rather than read from the LUT, so
% this compares against known values instead of the file against itself.
% RFPulseType sits in the printout block, which is where the PDF tool puts it
% and where a plain dcm2niix sidecar has nothing at all.
acq = struct( ...
    'ProtocolName',          'T2MAP_adapted_from_GRAPPATINI', ...
    'SequenceName',          'WIP_se2d16', ...
    'SoftwareVersions',      'syngo MR E11', ...
    'PulseSequenceDetails',  '%CustomerSeq%\accT2map_wip899c', ...
    'MagneticFieldStrength', 3, ...
    'SliceThickness',        2, ...
    'SiemensProtocolPrintout', struct('RFPulseType', 'Low SAR'));

rep = check_LUT_compatibility(G, acq);
check(numel(rep.Matched) == 7 && isempty(rep.Unverified), 'the shipped LUT matches all seven identity fields');

wrongs = { 'ProtocolName',          'T2MAP_something_else'
           'SequenceName',          'WIP_se2d32'
           'SoftwareVersions',      'syngo MR XA30'
           'PulseSequenceDetails',  '%CustomerSeq%\accT2map_wip900d'
           'MagneticFieldStrength', 1.5
           'SliceThickness',        3 };
for k = 1:size(wrongs,1)
    bad = acq;  bad.(wrongs{k,1}) = wrongs{k,2};
    check_error(@() check_LUT_compatibility(G, bad), 'check_LUT_compatibility:ProtocolMismatch', ...
        sprintf('a different %s is refused', wrongs{k,1}));
end
bad = acq;  bad.SiemensProtocolPrintout.RFPulseType = 'Normal';
check_error(@() check_LUT_compatibility(G, bad), 'check_LUT_compatibility:ProtocolMismatch', ...
    'a different RFPulseType is refused - the pulse is what the slice profile comes from');

% 2.89362 T and 3 T are the same magnet written two ways; 1.5 T is not, and is
% caught by the loop above.
near = acq;  near.MagneticFieldStrength = 2.89362;
accepted = true;
try
    repNear = check_LUT_compatibility(G, near);
catch
    accepted = false;      % a mismatch here is the failure this check is about
end
check(accepted && numel(repNear.Matched) == 7, '2.89362 T and 3 T count as the same field strength');

% A hand-written sidecar may put RFPulseType at the top level instead.
flat = rmfield(acq, 'SiemensProtocolPrintout');  flat.RFPulseType = 'Low SAR';
repFlat = check_LUT_compatibility_quietly(G, flat);
check(numel(repFlat.Matched) == 7, 'RFPulseType is found at the top level too');

check_error(@() check_LUT_compatibility(G, rmfield(acq,'ProtocolName')), 'check_LUT_compatibility:NoProtocolName', ...
    'an acquisition with no ProtocolName cannot be matched');
check_error(@() check_LUT_compatibility(struct('TrueT2',1), acq), 'check_LUT_compatibility:NoProtocol', ...
    'a struct that is not one of our LUTs is refused');
odd = acq;  odd.SliceThickness = [2 2 2];
check_error(@() check_LUT_compatibility(G, odd), 'check_LUT_compatibility:NotComparable', ...
    'a value that is neither text nor one number is refused rather than compared loosely');

% A thinner sidecar still matches on what it carries, and must say what it could not check.
thin = rmfield(acq, {'SequenceName','SoftwareVersions','SiemensProtocolPrintout'});
check_warning(@() check_LUT_compatibility(G, thin), 'check_LUT_compatibility:Unverified', ...
    'a sidecar missing identity fields warns rather than passing silently');
repThin = check_LUT_compatibility_quietly(G, thin);
check(isequal(repThin.Unverified, {'SequenceName','SoftwareVersions','RFPulseType'}), ...
    'and names the fields it could not check');

%% ------------------------------------------------------------------ %%
%  Worked example - the shipped table, end to end
%% ------------------------------------------------------------------ %%
fprintf('\nworked example (%s)\n', G.Protocol.ProtocolName);
fprintf('   LUT = jsondecode(fileread(fullfile(''LUT'',''%s'')));\n', shipped);
fprintf('   [correctedT2, info] = apply_T2_LUT(apparentT2, LUT, ''mask'', brainmask);\n\n');
demo = [25 60 90 130 420].';
[corr, dInfo] = apply_T2_LUT(demo, G);
fprintf('   %10s %12s %10s\n', 'apparent', 'corrected', 'change');
for k = 1:numel(demo)
    fprintf('   %10.1f %12.1f %10.1f\n', demo(k), corr(k), corr(k)-demo(k));
end
fprintf('   %d of %d voxels corrected, %d above the table (NaN)\n', ...
        dInfo.CorrectedVoxels, dInfo.ConsideredVoxels, dInfo.AboveRangeVoxels);
fprintf('   this table is valid for %s only - check Protocol against your sidecar\n', G.Protocol.ProtocolName);

fprintf('\nrelease_tests: all checks passed\n\n');

%% ------------------------------------------------------------------ %%
%  Helpers
%% ------------------------------------------------------------------ %%
function LUT = synthetic_LUT()
    % Same fields and shape as a real LUT, made up rather than simulated: a
    % monotone map with apparent > true and a slope slightly above 1. Used so
    % that the apply_T2_LUT checks do not depend on any particular table.
    trueT2 = 10:1:300;
    LUT = struct( ...
        'Description', 'Synthetic LUT for release_tests - not a real correction.', ...
        'T2Units',     'ms', ...
        'TrueT2',      trueT2, ...
        'ApparentT2',  trueT2 + 0.5 + 0.045*trueT2, ...
        'FitMethod',   'nls', ...
        'SkipEcho',    1, ...
        'T1_s',        1.33, ...
        'B1',          1);
end

function check(cond, name)
    if cond
        fprintf('   PASS  %s\n', name);
    else
        error('release_tests:Failed', 'FAIL  %s', name);
    end
end

function check_warning(f, expected_id, name)
    % Same idea as check_error, for a guard that warns instead of throwing.
    orig = warning('off', expected_id);
    restore = onCleanup(@() warning(orig)); %#ok<NASGU>
    lastwarn('', '');
    f();
    [~, got] = lastwarn;
    check(strcmp(got, expected_id), sprintf('%s [%s]', name, expected_id));
end

function r = check_LUT_compatibility_quietly(varargin)
    % The warning is the point of the check above; here only the report is.
    orig = warning('off', 'check_LUT_compatibility:Unverified');
    restore = onCleanup(@() warning(orig)); %#ok<NASGU>
    r = check_LUT_compatibility(varargin{:});
end

function check_error(f, expected_id, name)
    % Assert on the error identifier, never on a count of the word "error":
    % a count is satisfied by any guard firing, so the check stays green after
    % the specific guard it is named for is deleted.
    got = '(no error raised)';
    try
        f();
    catch ME
        got = ME.identifier;
    end
    check(strcmp(got, expected_id), sprintf('%s [%s]', name, expected_id));
end
