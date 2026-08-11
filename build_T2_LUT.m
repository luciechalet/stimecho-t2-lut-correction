function [LUT, S, TE] = build_T2_LUT(opt, args)
    arguments
        opt         (1,1) struct
        args.T2grid (1,:) double {mustBePositive} = 10:1:300            
        args.B1     (1,1) double {mustBeInRange(args.B1,0.1,2)} = 1     
        args.skipEcho = []
        args.method (1,:) char {mustBeMember(args.method,{'nls','loglinear'})} = 'nls'
        args.outdir (1,:) char = '' % build only, write nothing
    end
    % Build the T2 lookup table (LUT) for predefined slice profiles and acquisition parameters
    % Usage:
    %       [LUT, S, TE] = build_T2_LUT(opt)
    %       [LUT, S, TE] = build_T2_LUT(opt, 'T2grid', 10:1:300, 'B1', 1, 'skipEcho', [], 'method', 'nls') 
    % Inputs:
    %       opt: struct built from build_opt_from_protocol.m (Slice profile already generated with getRF.m)
    %       args.T2grid: range of true T2 values to simulate (in ms of the form: min:step:max, default: 10:1:300)
    %       args.B1: B1 field inhomogeneity, assumed constant and set to ideal 1, see assumptions in README for more information.
    %       args.skipEcho: Vector of echo indices to skip when applying the fit or integer for single echo (default: [])
    %       args.method: Fitting method, either 'nls' (non-linear least squares) or 'loglinear' (log-linear regression) (default: 'nls')
    %       args.outdir: directory to write the LUT JSON file (default: '' for build only, write nothing)
    % Outputs:
    %       LUT: Lookup table of apparent T2 values (ms) for each true T2 value in args.T2grid
    %       S: Simulated signal for each true T2 value in args.T2grid
    %       TE: Echo times for the simulated signal
    % Dependencies: StimFit toolbox
    require_stimfit(); % add StimFit to path if not already there
    %% Sanity check
    assert(isfield(opt,'RFe') && isfield(opt.RFe,'alpha'), 'build_T2_LUT:NoSliceProfiles', ...
        'opt carries no slice profiles - was it built by build_opt_from_protocol?');

    assert(numel(opt.RFe.alpha) == opt.Nz && numel(opt.RFr.alpha) == opt.Nz, 'build_T2_LUT:StaleSliceProfiles', ...
        'Slice profiles are stale or missing - re-run getRF(opt) after changing Dz/Nz.');
    
    assert(issorted(args.T2grid,'strictascend'), 'build_T2_LUT:UnsortedGrid', ...
        'T2grid must be strictly increasing - the apply step inverts it with interp1.');    

    assert(isfield(opt,'protocol') && isfield(opt.protocol,'ProtocolName'), 'build_T2_LUT:NoProtocol', ...
        'opt carries no protocol identity - rebuild it with build_opt_from_protocol.');


    %% Building T2 LUT
    T2grid = args.T2grid;
    TE     = (opt.esp:opt.esp:opt.esp*opt.etl).';
    S      = zeros(opt.etl, numel(T2grid));    

    for k = 1:numel(T2grid)
        S(:,k) = FSEsig(T2grid(k)*1e-3, args.B1, 0, opt);   % ms -> s; Gbg = 0 (see assumptions in README for more information)
    end

    %% Fit the monoexponential model to each simulated signal
    TEms = TE * 1000;              % build_T2_LUT returns TE in s; the LUT is in ms
    apparent_T2 = zeros(size(T2grid));
    for k = 1:numel(T2grid)
        apparent_T2(k) = fit_monoexp(S(:,k), TEms, 'method', args.method, 'skipEcho', args.skipEcho);
    end

    if ~issorted(apparent_T2,'strictascend')
        bad = find(diff(apparent_T2) <= 0, 1);
        error('build_T2_LUT:NotInvertible', ...
            ['Apparent T2 is not strictly increasing (true %g -> %.4f ms, then %g -> %.4f ms). ' ...
             'The apply step inverts this curve with interp1, which needs a one-to-one map.'], ...
            T2grid(bad), apparent_T2(bad), T2grid(bad+1), apparent_T2(bad+1));
    end

    %% Build LUT 
    T1_s = arrayfun(opt.T1, T2grid*1e-3); % T1 is a function of T2, convert ms -> s

    % Create LUT structure
    LUT = struct();
    LUT.Description = 'Apparent (monoexponential) vs True T2 from EPG simulation of the protocol.';
    LUT.Protocol = opt.protocol;
    LUT.T2Units = 'ms';
    LUT.TrueT2 = T2grid;
    LUT.ApparentT2 = apparent_T2;
    LUT.FitMethod = args.method;
    LUT.SkipEcho = parse_skipEcho(args.skipEcho);
    if all(T1_s == T1_s(1))
        LUT.T1_s = T1_s(1);
    else
        LUT.T1_s = T1_s;
    end
    LUT.B1 = args.B1;
    LUT.EchoSpacing_s = opt.esp;
    LUT.EchoTrainLength = opt.etl;
    LUT.SliceSamples = opt.Nz;
    LUT.RFwaveforms = {basename(opt.RFe.path), basename(opt.RFr.path)};
    LUT.Provenance = struct('stimecho_t2_lut_correction', git_info(fileparts(mfilename('fullpath'))),...
    'StimFit', git_info(fileparts(fileparts(which('FSEsig')))));
    LUT.Timestamp = char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));

    %% Write LUT to JSON file if requested
    if ~isempty(args.outdir)
        if ~isfolder(args.outdir)
            mkdir(args.outdir);
        end
        fpath = fullfile(args.outdir, lutfilename(LUT));
        fid = fopen(fpath, 'w');
        if fid < 0
            error('build_T2_LUT:CannotWrite', 'Could not open %s for writing.', fpath);
        end
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '%s', jsonencode(LUT, 'PrettyPrint', true));
        fprintf('LUT written to %s\n', fpath);
    end
end

function b = basename(p)
    [~, n, e] = fileparts(p);
    b = [n e];
end

function info = git_info(dir)
    % Get git information for the given directory, if it is a git repository.
    info = struct('Remote','', 'Branch','', 'Commit','', 'Dirty', false);
    if isempty(dir), return; end
    q = @(cmd) strtrim(system_out(sprintf('git -C "%s" %s', dir, cmd)));
    info.Commit = q('rev-parse --short HEAD');
    if isempty(info.Commit), return; end
    info.Branch = q('rev-parse --abbrev-ref HEAD');
    info.Dirty = ~isempty(q('status --porcelain --untracked-files=no'));

    % check remote
    remote_name =  q(sprintf('config --get branch.%s.remote', info.Branch));
    if isempty(remote_name), remote_name = 'origin'; end
    info.Remote = q(sprintf('remote get-url %s', remote_name));
end

function out = system_out(cmd)
    [st, out] = system(cmd);
    if st ~= 0, out = ''; end
end

function name = lutfilename(LUT)
    % Create LUT filename based on protocol name 
    if isscalar(LUT.T1_s)
        T1tag = sprintf('T1-%dms', round(LUT.T1_s*1000));
    else
        T1tag = 'T1-funcT2';
    end
    if isempty(LUT.SkipEcho)
        skiptag = 'skip-none';
    else
        skiptag = ['skip' char(strjoin(string(LUT.SkipEcho), '-'))];
    end
    B1tag = sprintf('B1-%03d', round(LUT.B1*100));
    protocoltag = regexprep(LUT.Protocol.ProtocolName, '[^\w-]', '_');
    name = sprintf('T2LUT_%s_%s_%s_%s_%s.json', protocoltag, T1tag, B1tag, skiptag, LUT.FitMethod);
end