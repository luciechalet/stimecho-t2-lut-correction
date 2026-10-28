% RELEASE_TESTS  Public test tier for stimecho-t2-lut-correction 
%
%   Run from the repository root:   matlab -batch release_tests
%
%   Covers everything that runs without proprietary data: apply_T2_LUT,
%   fit_monoexp, parse_skipEcho, and the LUT files shipped under LUT/.

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
% monoexponential curve
[t2_nls, s0_nls] = fit_monoexp(S_pure, TE_ms, 'method', 'nls');
[t2_log, s0_log] = fit_monoexp(S_pure, TE_ms, 'method', 'loglinear');
check(abs(t2_nls - T2_true) < 1e-6,  sprintf('nls recovers a pure exponential (%.10f ms)', t2_nls));
check(abs(t2_log - T2_true) < 1e-9,  sprintf('loglinear recovers it in closed form (%.12f ms)', t2_log));
check(abs(s0_nls - S0_true) < 1e-5 && abs(s0_log - S0_true) < 1e-8, 'both recover S0');

% Unit-agnostic: T2 comes back in whatever units TE was given in.
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
% and not the contents of any particular LUT.
LUT = synthetic_LUT();
appAxis = LUT.ApparentT2(:);
trueAxis= LUT.TrueT2(:);

% Round trip: feeding the table its own apparent axis must return its own true
% axis.
back = apply_T2_LUT(appAxis, LUT);
check(max(abs(back - trueAxis)) < 1e-9, sprintf('round trip through the table (max %.2e ms)', max(abs(back - trueAxis))));
check(all(back < appAxis), 'the correction moves T2 DOWN (apparent > true)');

% Out of range is NaN at both ends: no clipping, no extrapolation.
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
%  apply_T2_LUT with T1 and B1 axes - base MATLAB only
%% ------------------------------------------------------------------ %%
fprintf('\napply_T2_LUT, T1 and B1 axes\n');

L3 = jsondecode(jsonencode(synthetic_LUT3()));   % through JSON, like a real file
t3 = L3.TrueT2(:);
a3 = synthetic_LUT3().ApparentT2(:,2,2);         % the curve at T1 = 1.3 s, B1 = 0.9
on = @(v) v*ones(size(a3));
inner = t3 >= 20 & t3 <= 290;                    % the ends can be NaN by design (edge of neighbouring curves)

% Round trip on a grid point: the correction must return the true axis.
c3 = apply_T2_LUT(a3, L3, 'T1map', on(1.3), 'B1map', on(0.9));
check(max(abs(c3(inner) - t3(inner))) < 1e-6 && ~any(isnan(c3(inner))), ...
    sprintf('round trip on a (T1, B1) grid point (max %.2e ms)', max(abs(c3(inner) - t3(inner)))));

% The maps must actually be used: same apparent T2, different B1, different answer.
cB = apply_T2_LUT(a3, L3, 'T1map', on(1.3), 'B1map', on(1.0));
check(min(abs(cB(inner) - c3(inner))) > 0.1, 'a different B1 gives a different correction');

% Outside the grid is NaN, and counted.
[cO, iO] = apply_T2_LUT(a3, L3, 'T1map', on(1.3), 'B1map', on(0.5));
check(all(isnan(cO)) && iO.OutsideB1GridVoxels == numel(a3), 'B1 outside the grid gives NaN and is counted');

% A single-value axis is dropped, not interpolated over.
L1 = L3;  L1.T1_s = 1.3;  L1.ApparentT2 = L3.ApparentT2(:,2,:);
c1 = apply_T2_LUT(a3, L1, 'B1map', on(0.9));
check(max(abs(c1(inner) - c3(inner))) < 1e-9, 'a LUT with one T1 value needs only the B1 map');

% T1 given as a function of T2 stores one T1 per T2 in T1_s: a record, not an axis.
Lf = synthetic_LUT();  Lf.T1_s = 0.5*Lf.TrueT2/1000 + 1.5;
check(isequal(apply_T2_LUT(80, Lf), apply_T2_LUT(80, synthetic_LUT())), 'a per-T2 T1_s is read as a record, not an axis');
check_error(@() apply_T2_LUT(80, Lf, 'T1map', 1.3), 'apply_T2_LUT:UnusedMap', 'and a T1 map is refused for it');

% Guards.
check_error(@() apply_T2_LUT(a3, L3, 'B1map', on(0.9)), 'apply_T2_LUT:MissingMap', 'a T1 axis without a T1 map is refused');
check_error(@() apply_T2_LUT(a3, L1, 'T1map', on(1.3), 'B1map', on(0.9)), 'apply_T2_LUT:UnusedMap', 'a map the LUT has no axis for is refused, not ignored');
check_error(@() apply_T2_LUT(80, synthetic_LUT(), 'B1map', 1), 'apply_T2_LUT:UnusedMap', 'maps passed to a 1-D LUT are refused');
check_error(@() apply_T2_LUT(a3, L3, 'T1map', on(1.3), 'B1map', 0.9), 'apply_T2_LUT:MapSizeMismatch', 'a map must be on the T2 grid');
check_error(@() apply_T2_LUT(a3, L3, 'T1map', on(1300), 'B1map', on(0.9)), 'apply_T2_LUT:T1Units', 'a T1 map in ms is rejected');
check_error(@() apply_T2_LUT(a3, L3, 'T1map', on(1.3), 'B1map', on(90)), 'apply_T2_LUT:B1Units', 'a B1 map in percent is rejected');

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

% Probe values for the table the repository ships. 
shipped = 'T2LUT_T2MAP_adapted_from_GRAPPATINI_T1-1330ms_B1-100_skip1_nls.json';
check(isfile(fullfile(here,'LUT',shipped)), 'the GRAPPATINI table is the one that ships');
G = jsondecode(fileread(fullfile(here,'LUT',shipped)));

check(strcmp(G.FitMethod,'nls') && isequal(G.SkipEcho,1) && G.B1 == 1 && abs(G.T1_s - 1.33) < 1e-12, ...
      'built with nls, echo 1 skipped, B1 = 1, T1 = 1.33 s');
check(isequal(G.TrueT2(:).', 10:1:300), 'grid is 10:1:300 ms, 291 points');

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

acq = struct( ...
    'ProtocolName',          'T2MAP_adapted_from_GRAPPATINI', ...
    'SequenceName',          'WIP_se2d16', ...
    'SoftwareVersions',      'syngo MR E11', ...
    'PulseSequenceDetails',  '%CustomerSeq%\accT2map_wip899c', ...
    'MagneticFieldStrength', 3, ...
    'SliceThickness',        2, ...
    'EchoTrainLength',       16, ...
    'EchoTime',              (1:16).' * 0.0111, ...
    'SiemensProtocolPrintout', struct('RFPulseType', 'Low SAR'));

rep = check_LUT_compatibility(G, acq);
check(numel(rep.Matched) == 9 && isempty(rep.Unverified), ...
    'the shipped LUT matches all seven identity fields and the echo train');

% The echo times are the run-01 / run-02 case: everything else identical,
% only the spacing differs - 13.3 ms against the table's 11.1 ms.
wrongs = { 'ProtocolName',          'T2MAP_something_else'
           'SequenceName',          'WIP_se2d32'
           'SoftwareVersions',      'syngo MR XA30'
           'PulseSequenceDetails',  '%CustomerSeq%\accT2map_wip900d'
           'MagneticFieldStrength', 1.5
           'SliceThickness',        3
           'EchoTrainLength',       12
           'EchoTime',              (1:16).' * 0.0133 };
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
check(accepted && numel(repNear.Matched) == 9, '2.89362 T and 3 T count as the same field strength');

% A hand-written sidecar may put RFPulseType at the top level instead.
flat = rmfield(acq, 'SiemensProtocolPrintout');  flat.RFPulseType = 'Low SAR';
repFlat = check_LUT_compatibility_quietly(G, flat);
check(numel(repFlat.Matched) == 9, 'RFPulseType is found at the top level too');

% The DICOM route records the pulse type as the scanner's code (RFPulseTypeCode),
% never as the label. A code is compared with a code only.
Gc = G;  Gc.Protocol = rmfield(G.Protocol, 'RFPulseType');  Gc.Protocol.RFPulseTypeCode = 4;
coded = rmfield(acq, 'SiemensProtocolPrintout');  coded.RFPulseTypeCode = 4;
repCode = check_LUT_compatibility_quietly(Gc, coded);
check(numel(repCode.Matched) == 9 && isempty(repCode.Unverified), ...
    'a LUT with RFPulseTypeCode matches an acquisition with the same code');
coded.RFPulseTypeCode = 2;
check_error(@() check_LUT_compatibility(Gc, coded), 'check_LUT_compatibility:ProtocolMismatch', ...
    'a different RFPulseTypeCode is refused');
repMixed = check_LUT_compatibility_quietly(Gc, acq);    % LUT has the code, acquisition the label
check(isequal(repMixed.Unverified, {'RFPulseType','RFPulseTypeCode'}), ...
    'a label is never compared with a code - both are reported unverified');

% Echo train. A dcm2niix sidecar of a GRAPPATINI map has no EchoTime at all, and a
% single-echo image has one TE: neither says anything about the train.
noTE = rmfield(acq, 'EchoTime');
repNoTE = check_LUT_compatibility_quietly(G, noTE);
check(isequal(repNoTE.Unverified, {'EchoTime'}), 'without echo times, the echo train is reported unverified');
oneTE = acq;  oneTE.EchoTime = 0.08;
repOneTE = check_LUT_compatibility_quietly(G, oneTE);
check(isequal(repOneTE.Unverified, {'EchoTime'}), 'a single TE is not mistaken for a train');
inMs = acq;  inMs.EchoTime = acq.EchoTime * 1000;
repMs = check_LUT_compatibility_quietly(G, inMs);
check(numel(repMs.Matched) == 9, 'echo times written in ms are read as ms');
shifted = acq;  shifted.EchoTime = acq.EchoTime + 0.002;
check_error(@() check_LUT_compatibility(G, shifted), 'check_LUT_compatibility:ProtocolMismatch', ...
    'the right spacing with the first echo offset is refused - the simulation puts it at one spacing');

% The acquisition may be given as a JSON file instead of a struct.
jf = [tempname '.json'];
fid = fopen(jf, 'w');  fprintf(fid, '%s', jsonencode(acq));  fclose(fid);
repFile = check_LUT_compatibility(G, jf);
delete(jf);
check(numel(repFile.Matched) == 9, 'a JSON file path is read like the struct it contains');
check_error(@() check_LUT_compatibility(G, 'no_such_file.json'), 'check_LUT_compatibility:BadAcquisition', ...
    'a path that is neither a JSON file nor a folder is refused');
check_error(@() check_LUT_compatibility(G, acq, 'SeriesNumber', 7), 'check_LUT_compatibility:SeriesWithoutDICOM', ...
    'SeriesNumber with a struct is refused rather than ignored');

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
check(isequal(repThin.Unverified, {'SoftwareVersions','SequenceName','RFPulseType','RFPulseTypeCode'}), ...
    'and names the fields it could not check');

% syngo XA exports no SequenceName, so a LUT built from XA data has none: that LUT
% must still be checkable, with SequenceName reported unverified.
Gxa = G;  Gxa.Protocol = rmfield(G.Protocol, 'SequenceName');
repXa = check_LUT_compatibility_quietly(Gxa, acq);
check(numel(repXa.Matched) == 8 && isequal(repXa.Unverified, {'SequenceName'}), ...
    'a LUT without SequenceName (syngo XA) is checked on everything else');

%% ------------------------------------------------------------------ %%
%  Siemens DICOM protocol reader - base MATLAB only
%% ------------------------------------------------------------------ %%
fprintf('\nread_siemens_ascconv / siemens_protocol_facts\n');

% A synthetic protocol block between binary bytes, as it sits inside a DICOM file.
asc = {
    'tSequenceFileName	 = 	""%CustomerSeq%\wip_899_AccT2_se_mc""'
    'lContrasts	 = 	3'
    'alTR[0]	 = 	5040000'
    'alTE[0]	 = 	10100'
    'alTE[1]	 = 	20200'
    'alTE[2]	 = 	30300'
    'alTE[3]	 = 	1700'
    'alTE[4]	 = 	1800'
    'adFlipAngleDegree[0]	 = 	180.0'
    'sSliceArray.asSlice[0].dThickness	 = 	4.0'
    'sTXSPEC.ucRFPulseType	 = 	4'
    };
fixture = [tempname '.dcm'];
fid = fopen(fixture, 'w');
fwrite(fid, uint8([0 1 2 255]));
fwrite(fid, uint8(['### ASCCONV BEGIN object=MrProtDataImpl@MrProtocolData ###' newline ...
                   strjoin(asc, newline) newline '### ASCCONV END ###']));
fwrite(fid, uint8([0 0 7]));
fclose(fid);
[P, found] = read_siemens_ascconv(fixture);
delete(fixture);

check(found && strcmp(P('alTE[1]'), '20200'), 'the block is found among binary bytes and read as text');
check(strcmp(P('sSliceArray.asSlice[0].dThickness'), '4.0'), 'dotted and indexed names are kept whole');

facts = siemens_protocol_facts(P);
check(isequal(facts.EchoTime, [0.0101; 0.0202; 0.0303]), ...
    'only the first lContrasts alTE values are echoes, converted us -> s');
check(facts.EchoTrainLength == 3 && facts.RepetitionTime == 5.04, 'echo count from lContrasts, TR us -> s');
check(facts.FlipAngle == 180 && facts.SliceThickness == 4, 'flip angle and slice thickness');
check(strcmp(facts.PulseSequenceDetails, '%CustomerSeq%\wip_899_AccT2_se_mc'), ...
    'the sequence path loses its doubled quotes');
check(facts.RFPulseTypeCode == 4, 'the RF pulse type is kept as the scanner code');

noEcho = containers.Map(P.keys, P.values);  remove(noEcho, 'lContrasts');
check_error(@() siemens_protocol_facts(noEcho), 'siemens_protocol_facts:MissingParameter', ...
    'a protocol without lContrasts is refused, not read as "all alTE values"');
past = containers.Map(P.keys, P.values);  past('lContrasts') = '5';
check_error(@() siemens_protocol_facts(past), 'siemens_protocol_facts:EchoesNotIncreasing', ...
    'reading past the real echoes into the template values is caught');

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

function LUT = synthetic_LUT3()
    % Same fields and shape as a LUT with T1 and B1 axes, made up rather than
    % simulated. Each (T1, B1) curve is a monotone line with apparent > true, and
    % its slope changes with T1 and B1, so reading the wrong curve gives a
    % different answer.
    trueT2 = (10:1:300).';
    T1 = [0.8 1.3 4.3];
    B1 = [0.8 0.9 1.0];
    App = zeros(numel(trueT2), numel(T1), numel(B1));
    for j = 1:numel(T1)
        for b = 1:numel(B1)
            App(:,j,b) = trueT2 + 0.5 + (0.03 + 0.01*T1(j) + 0.1*(1 - B1(b)))*trueT2;
        end
    end
    LUT = struct( ...
        'Description',    'Synthetic 3-D LUT for release_tests - not a real correction.', ...
        'T2Units',        'ms', ...
        'TrueT2',         trueT2.', ...
        'ApparentT2',     App, ...
        'FitMethod',      'nls', ...
        'SkipEcho',       1, ...
        'T1_s',           T1, ...
        'B1',             B1, ...
        'ApparentT2Dims', {{'TrueT2', 'T1_s', 'B1'}});   % {{ }}: one cell as the field value
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
    % Assert on the error identifier
    got = '(no error raised)';
    try
        f();
    catch ME
        got = ME.identifier;
    end
    check(strcmp(got, expected_id), sprintf('%s [%s]', name, expected_id));
end
