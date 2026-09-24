function [LUT, S, TE] = build_T2_LUT(opt, args)
    arguments
        opt         (1,1) struct
        args.T2grid (1,:) double {mustBePositive} = 10:1:300
        args.T1grid (1,:) double {mustBePositive} = []   % seconds, [] = use the T1 already in opt (1-D LUT as before)
        args.B1     (1,:) double {mustBeInRange(args.B1,0.1,2)} = 1   % scalar, or a grid for a B1 axis
        args.skipEcho = []
        args.method (1,:) char {mustBeMember(args.method,{'nls','loglinear'})} = 'nls'
        args.outdir (1,:) char = '' % build only, write nothing
        args.parallel (1,1) logical = false   % true: spread simulations over a parallel pool (Parallel Computing Toolbox)
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

    assert(issorted(args.T1grid,'strictascend') && issorted(args.B1,'strictascend'), 'build_T2_LUT:UnsortedGrid', ...
        'T1grid and B1 must be strictly increasing - the apply step interpolates along them.');

    assert(all(args.T1grid < 20), 'build_T2_LUT:T1Units', ...
        'T1grid looks like milliseconds; it must be in seconds.');

    assert(isfield(opt,'protocol') && isfield(opt.protocol,'ProtocolName'), 'build_T2_LUT:NoProtocol', ...
        'opt carries no protocol identity - rebuild it with build_opt_from_protocol.');


    %% Building T2 LUT: simulate and fit for every (T2, T1, B1) combination
    T2grid   = args.T2grid;
    T1grid   = args.T1grid;
    B1grid   = args.B1;
    method   = args.method;       % plain variables rather than args.x: parfor sends them to each worker once
    skipEcho = args.skipEcho;
    TE       = (opt.esp:opt.esp:opt.esp*opt.etl).';
    TEms     = TE * 1000;         % the LUT is in ms
    nT2  = numel(T2grid);
    nT1  = max(1, numel(T1grid)); % 1 when T1 comes from opt
    nB1  = numel(B1grid);
    nSim = nT2 * nT1 * nB1;

    nWorkers = 0;                 % 0 = ordinary serial loop, no pool
    if args.parallel
        assert(license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel')), ...
            'build_T2_LUT:NoParallel', 'parallel = true needs the Parallel Computing Toolbox.');
        pool = gcp;               % reuse the running pool, or start one
        nWorkers = pool.NumWorkers;
    end

    % One flat loop over all combinations, so parfor can hand them out freely
    S_all   = zeros(opt.etl, nSim);
    app_all = zeros(1, nSim);
    parfor (i = 1:nSim, nWorkers)
        [k, j, b] = ind2sub([nT2 nT1 nB1], i);   % flat index -> (T2, T1, B1) indices
        optj = opt;
        if ~isempty(T1grid)
            optj.T1 = @(T2) T1grid(j);           % same pattern as build_opt_from_protocol
        end
        s = FSEsig(T2grid(k)*1e-3, B1grid(b), 0, optj);   % ms -> s; Gbg = 0 (see README)
        S_all(:,i) = s;
        app_all(i) = fit_monoexp(s, TEms, 'method', method, 'skipEcho', skipEcho);
    end
    S           = reshape(S_all,   opt.etl, nT2, nT1, nB1);
    apparent_T2 = reshape(app_all, nT2, nT1, nB1);

    % Each (T1, B1) curve must be invertible on its own
    for j = 1:nT1
        for b = 1:nB1
            col = apparent_T2(:,j,b);
            if ~issorted(col, 'strictascend')
                bad = find(diff(col) <= 0, 1);
                if isempty(T1grid), T1val = opt.T1(T2grid(bad)*1e-3); else, T1val = T1grid(j); end
                error('build_T2_LUT:NotInvertible', ...
                    ['Apparent T2 is not strictly increasing at T1 = %g s, B1 = %g ' ...
                     '(true %g -> %.4f ms, then %g -> %.4f ms).'], ...
                    T1val, B1grid(b), T2grid(bad), col(bad), T2grid(bad+1), col(bad+1));
            end
        end
    end

    %% Create LUT structure
    LUT = struct();
    LUT.Description = 'Apparent (monoexponential) vs True T2 from EPG simulation of the protocol.';
    LUT.Protocol = opt.protocol;
    LUT.T2Units = 'ms';
    LUT.TrueT2 = T2grid;
    LUT.ApparentT2 = apparent_T2;
    LUT.FitMethod = args.method;
    LUT.SkipEcho = parse_skipEcho(args.skipEcho);
    if isempty(args.T1grid)
        T1_s = arrayfun(opt.T1, T2grid*1e-3);    % T1 from opt: one value, or one per T2
        if all(T1_s == T1_s(1))
            LUT.T1_s = T1_s(1);
        else
            LUT.T1_s = T1_s;
        end
    else
        LUT.T1_s = args.T1grid;                  % T1 axis
    end
    % Axes of ApparentT2 that are grids. ApparentT2 is always ordered (TrueT2, T1, B1);
    % a T1 or B1 not listed here has a length-1 dimension, and T1_s then records the T1 used.
    LUT.ApparentT2Dims = {'TrueT2'};
    if ~isempty(args.T1grid), LUT.ApparentT2Dims{end+1} = 'T1_s'; end
    if numel(args.B1) > 1,    LUT.ApparentT2Dims{end+1} = 'B1';   end
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
    if any(strcmp('T1_s', LUT.ApparentT2Dims))
        T1tag = sprintf('T1-%d-%dms', round(LUT.T1_s([1 end])*1000));
    elseif isscalar(LUT.T1_s)
        T1tag = sprintf('T1-%dms', round(LUT.T1_s*1000));
    else
        T1tag = 'T1-funcT2';
    end
    if isempty(LUT.SkipEcho)
        skiptag = 'skip-none';
    else
        skiptag = ['skip' char(strjoin(string(LUT.SkipEcho), '-'))];
    end
    if isscalar(LUT.B1)
        B1tag = sprintf('B1-%03d', round(LUT.B1*100));
    else
        B1tag = sprintf('B1-%03d-%03d', round(LUT.B1([1 end])*100));
    end
    protocoltag = regexprep(LUT.Protocol.ProtocolName, '[^\w-]', '_');
    name = sprintf('T2LUT_%s_%s_%s_%s_%s.json', protocoltag, T1tag, B1tag, skiptag, LUT.FitMethod);
end