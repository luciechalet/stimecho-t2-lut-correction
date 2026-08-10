function idx = parse_skipEcho(spec)
    % This function parses the skipEcho string and returns the indices of the echoes to skip.
    % The skipEcho string can be a comma-separated list of indices (e.g "1,3,5").
    % The function returns a vector of indices to skip

    %   idx = parse_skipEcho("")    -> []        (skip nothing)
    %   idx = parse_skipEcho("1")   -> 1
    %   idx = parse_skipEcho("1,3") -> [1 3]
    %   idx = parse_skipEcho(1)     -> 1       

    if isnumeric(spec)
        idx = double(spec(:)).';
    else
        spec = strtrim(string(spec));
        if spec == ""
            idx = [];
        else
            idx = str2double(split(spec, ',')).';
        end
    end

    if isempty(idx)
        idx = [];
    end

    if any(~isfinite(idx)) || any(idx < 1) || any(mod(idx,1) ~= 0)
        error('parse_skipEcho:BadSpec', ...
            ['Invalid skipEcho. Expected a non negative integer or a comma-separated list of non-negative integers (e.g. "1,3,5").' ...
            ' or an empty string; got "%s", parsed as "%s".'], string(spec), mat2str(idx));
    end
    end