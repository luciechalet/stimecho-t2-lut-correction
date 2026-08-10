function [T2, S0] = fit_monoexp(S, TE, args)
    arguments
        S           (:,1) double % Signal 
        TE          (:,1) double {mustBeNonnegative} % Echo times
        args.method (1,:) char {mustBeMember(args.method,{'nls','loglinear'})} = 'nls'  % Fitting method
        args.skipEcho                                                          = [] % Echoes to skip when applying fit
    end
    % Fit a monoexponential |S| = A*exp(-TE/T2) 
    %  [T2, S0] = fit_monoexp(S, TE, 'method', 'nls', 'skipEcho', []); For non least squares fitting (default) without skipping any echoes (default)
    %  [T2, S0] = fit_monoexp(S, TE, 'method', 'nls', 'skipEcho', [1]); For non least squares fitting (default) while skipping the first echo (GRAPPATINI)
    %  [T2, S0] = fit_monoexp(S, TE, 'method', 'loglinear', 'skipEcho', []); For log-linear fitting without skipping any echoes (default)

    % Inputs:
    %  S: Signal vector (amplitude) of size (nEchoes,1)
    %  TE: Echo times vector of size (nEchoes,1)
    %  args.method: Fitting method, either 'nls' (non-linear least squares) or 'loglinear' (log-linear regression) (default: 'nls')
    %  args.skipEcho: Vector of echo indices to skip when applying the fit (default: [])

    % Outputs:
    %  T2: Fitted T2 value, in the same units as TE
    %  S0: Fitted S0 value (amplitude at TE=0)

    %% Select echoes to fit
    if numel(S) ~= numel(TE)
        error('fit_monoexp:SizeMismatch', ...
            'S has %d echoes but TE has %d.', numel(S), numel(TE));
    end

    skip = parse_skipEcho(args.skipEcho);
    if any(skip > numel(TE))
        error('fit_monoexp:SkipOutOfRange', ...
            ['skipEcho asks for echo %d but the train has only %d echoes. ' ...
             'setdiff would drop it silently, so it is caught here.'], ...
            max(skip), numel(TE));
    end
    keep = setdiff(1:numel(TE), skip);      % setdiff(1:n, []) is 1:n
    if numel(keep) < 3
        error('fit_monoexp:TooFewEchoes', ...
            ['Only %d echo(es) left after skipEcho - a two-parameter fit needs ' ...
             'at least 3 to be overdetermined.'], numel(keep));
    end

    t = TE(keep);
    y = abs(S(keep));

    %% Fit the monoexponential model
    % initial T2 guess using log linear fit
    p = polyfit(t, log(y), 1);
    switch args.method
        case 'loglinear'
            % Log-linear regression: log(S) = log(S0) - TE/T2
            T2 = -1/p(1);
            S0 = exp(p(2));
        case 'nls'
            % Non-linear least squares fitting: S = S0 * exp(-TE/T2)
            % Initial guess for T2 using the first log linear estimate
            T2_init = -1/p(1);
            
            % Set fitting options
            fopt = optimset('lsqnonlin');
            fopt.Display = 'off';
            fopt.TolX    = 1e-8;
            fopt.TolFun  = 1e-9;

            % Applying unbounded non-linear least squares fitting using stimfit's lsqnonlin
            [B, ~, ~, exitflag] = lsqnonlin(@(x) x(2)*exp(-t/x(1)) - y, [T2_init, max(y)*1.5], [], [], fopt); 

            % Considering a noiseless simulation, failed convergence is unlikely and throws an error
            if exitflag <= 0
                error('fit_monoexp:NoConvergence', ...
                    'lsqnonlin did not converge (exitflag %d, init %g).', ...
                    exitflag, T2_init);
            end
            T2 = B(1);
            S0 = B(2);
    end
end

