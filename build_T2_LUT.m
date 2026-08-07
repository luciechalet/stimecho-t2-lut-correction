function [S, TE, T2grid] = build_T2_LUT(opt, args)
    arguments
        opt         (1,1) struct
        args.T2grid (1,:) double {mustBePositive} = 20:1:200            % ms, D2
        args.B1     (1,1) double {mustBeInRange(args.B1,0.1,2)} = 1     % D5
    end
    % Build the T2 lookup table (LUT) for predifined slice profiles and acquisition parameters
    %   [S, TE, T2grid] = build_T2_LUT(opt)
    %   [S, TE, T2grid] = build_T2_LUT(opt, 'T2grid', 20:1:200, 'B1', 1) 

    % Dependencies: StimFit toolbox
    require_stimfit(); % add StimFit to path if not already there

    % Inputs:
    %  opt: struct built from build_opt_from_protocol.m (Slice profile already generated with getRF.m)
    %  args.T2grid: range of true T2 values to simulate (in ms of the form: min:step:max, default: 20:1:200)
    %  args.B1: B1 field inhomogeneity, assumed constant and set to ideal 1, see assumptions in README for more information.

    % Outputs:
    %  TODO: update later

    %% Sanity check
    assert(isfield(opt,'RFe') && isfield(opt.RFe,'alpha'), ...
        'opt carries no slice profiles - was it built by build_opt_from_protocol?');

    assert(numel(opt.RFe.alpha) == opt.Nz && numel(opt.RFr.alpha) == opt.Nz, ...
        'Slice profiles are stale or missing - re-run getRF(opt) after changing Dz/Nz.');
    
    assert(issorted(args.T2grid,'strictascend'), ...
        'T2grid must be strictly increasing - the apply step inverts it with interp1.');    

    %% Building T2 LUT
    T2grid = args.T2grid;
    TE     = (opt.esp:opt.esp:opt.esp*opt.etl).';
    S      = zeros(opt.etl, numel(T2grid));    

    for k = 1:numel(T2grid)
        S(:,k) = FSEsig(T2grid(k)*1e-3, args.B1, 0, opt);   % ms -> s; Gbg = 0 (see assumptions in README for more information)
    end
end
