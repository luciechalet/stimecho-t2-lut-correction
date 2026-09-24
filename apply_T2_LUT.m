function [correctedT2, info] = apply_T2_LUT(apparentT2, LUT, args)
    arguments 
        apparentT2 {mustBeNumeric}
        LUT (1,1) struct
        args.mask {mustBeNumericOrLogical} = []
        args.T1map {mustBeNumeric} = []   % seconds, same size as apparentT2, required if the LUT has a T1 axis
        args.B1map {mustBeNumeric} = []   % flip-angle ratio (1 = nominal), same size, required if the LUT has a B1 axis
    end
    % Apply the T2 LUT to correct stimulated echoes in apparent T2 values.
    % Usage:
    %       [correctedT2, info] = apply_T2_LUT(apparentT2, LUT);
    %       [correctedT2, info] = apply_T2_LUT(apparentT2, LUT, 'mask', ROI_mask);
    % Inputs:
    %       apparentT2: Apparent T2 values array (must be in ms) 
    %       LUT: Look up table structure built with build_T2_LUT
    %       args.mask: Optional logical or numeric mask to apply the correction only to selected voxels (voxels outside the mask will be returned as NaN)
    % Outputs:
    %       correctedT2: Corrected T2 values array (in ms) same dimensions as apparentT2
    %       info: Structure containing information about the LUT and the correction process

    %% Sanity checks: LUT structure
    required_fields = {'TrueT2', 'ApparentT2', 'T2Units', 'B1'};
    missing_fields = required_fields(~isfield(LUT, required_fields));
    assert(isempty(missing_fields), 'apply_T2_LUT:MissingField', ...
        'LUT is missing field(s): %s. Expected a struct from jsondecode of a build_T2_LUT JSON.', ...
        strjoin(missing_fields, ', '));
    assert(strcmpi(LUT.T2Units, 'ms'), 'apply_T2_LUT:WrongUnits', ...
        'LUT.T2Units is %s, but apply_T2_LUT expects LUT in ms.', LUT.T2Units);

    x = double(apparentT2); % ensure double for interpolation
    T1isAxis = isfield(LUT, 'ApparentT2Dims') && any(strcmp('T1_s', LUT.ApparentT2Dims));
    hasAxes  = T1isAxis || numel(LUT.B1) > 1;

    if ~hasAxes
        %% 1-D LUT: unchanged behaviour
        assert(isempty(args.T1map) && isempty(args.B1map), 'apply_T2_LUT:UnusedMap', ...
            'This LUT has no T1 or B1 axis, so the maps would be ignored.');
        apparentAxis = LUT.ApparentT2(:);
        trueAxis     = LUT.TrueT2(:);
        assert(numel(apparentAxis) == numel(trueAxis), 'apply_T2_LUT:LengthMismatch', ...
            'LUT.ApparentT2 (%d) and LUT.TrueT2 (%d) differ in length.', ...
            numel(apparentAxis), numel(trueAxis));
        assert(all(diff(apparentAxis) > 0), 'apply_T2_LUT:NotInvertible', ...
            'LUT.ApparentT2 is not strictly increasing, so the table cannot be inverted.');
        correctedT2 = interp1(apparentAxis, trueAxis, x, 'pchip', NaN); % NaN for out-of-bounds values, pchip for shape preserving cubic interpolation
    else
        %% LUT with a T1 and/or B1 axis: per-voxel maps
        [correctedT2, apparentAxis] = apply_multidim(x, LUT, args.T1map, args.B1map);
    end

    %% Apply mask if provided
    if ~isempty(args.mask)
        assert(isequal(size(args.mask), size(apparentT2)), 'apply_T2_LUT:MaskSizeMismatch', ...
            'Mask size (%s) does not match apparentT2 size (%s).', ...
            mat2str(size(args.mask)), mat2str(size(apparentT2)));
        keep = logical(args.mask);
    else
        keep = true(size(apparentT2));
    end
    correctedT2(~keep) = NaN; % set voxels outside the mask to NaN

    %% Prepare info structure
    consideredVoxels = keep & isfinite(x);

    info = struct( ...
        'ConsideredVoxels', nnz(consideredVoxels), ...
        'CorrectedVoxels', nnz(consideredVoxels & isfinite(correctedT2)), ...
        'BelowRangeVoxels', nnz(consideredVoxels & (x < apparentAxis(1))), ...
        'AboveRangeVoxels', nnz(consideredVoxels & (x > apparentAxis(end))), ...
        'ApparentT2Range', [apparentAxis(1), apparentAxis(end)], ...
        'T2Units', 'ms');
    
    if ~isempty(args.T1map)
        info.OutsideT1GridVoxels = nnz(consideredVoxels & ...
            (args.T1map < LUT.T1_s(1) | args.T1map > LUT.T1_s(end)));
    end
    if ~isempty(args.B1map)
        info.OutsideB1GridVoxels = nnz(consideredVoxels & ...
            (args.B1map < LUT.B1(1) | args.B1map > LUT.B1(end)));
    end

end

function [correctedT2, appAxis] = apply_multidim(x, LUT, T1map, B1map)
    % Correct with a LUT that has a T1 and/or B1 axis, using per-voxel maps.
    trueAxis = LUT.TrueT2(:);
    B1axis   = LUT.B1(:);
    if isfield(LUT, 'ApparentT2Dims') && any(strcmp('T1_s', LUT.ApparentT2Dims))
        T1axis = LUT.T1_s(:);
    else
        T1axis = NaN;
    end
    useT1 = numel(T1axis) > 1;   % an axis with a single value is not interpolated over
    useB1 = numel(B1axis) > 1;

    % Each map is required exactly when the LUT has an axis for it
    check_map(T1map, useT1, 'T1', size(x));
    check_map(B1map, useB1, 'B1', size(x));
    if useT1
        medT1 = median(T1map(T1map > 0 & isfinite(T1map)));
        assert(medT1 < 20, 'apply_T2_LUT:T1Units', ...
            'T1map median is %g; it looks like milliseconds, it must be in seconds.', medT1);
    end
    if useB1
        medB1 = median(B1map(B1map > 0 & isfinite(B1map)));
        assert(medB1 > 0.3 && medB1 < 2, 'apply_T2_LUT:B1Units', ...
            'B1map median is %g; expected a flip-angle ratio near 1 (not percent, not degrees x10).', medB1);
    end

    % Restore the (TrueT2, T1, B1) shape in case the JSON round trip changed it
    App = reshape(LUT.ApparentT2, numel(trueAxis), numel(T1axis), numel(B1axis));

    % Invert every (T1, B1) curve onto one shared apparent-T2 axis
    step    = 0.05;   % ms
    appAxis = unique([(min(App(:)) : step : max(App(:))).'; max(App(:))]);   % make sure the axis reaches the top value
    Inv = NaN(numel(appAxis), numel(T1axis), numel(B1axis));
    for j = 1:numel(T1axis)
        for b = 1:numel(B1axis)
            col = App(:,j,b);
            assert(all(diff(col) > 0), 'apply_T2_LUT:NotInvertible', ...
                'LUT curve at T1 index %d, B1 index %d is not strictly increasing.', j, b);
            Inv(:,j,b) = interp1(col, trueAxis, appAxis, 'pchip', NaN);
        end
    end

    % One lookup per voxel at (apparent T2, T1, B1); single-value axes are dropped
    grids  = {appAxis};
    points = {x(:)};
    if useT1, grids{end+1} = T1axis; points{end+1} = double(T1map(:)); end
    if useB1, grids{end+1} = B1axis; points{end+1} = double(B1map(:)); end
    F = griddedInterpolant(grids, squeeze(Inv), 'linear', 'none');   % 'none' = NaN outside the grids
    correctedT2 = reshape(F(points{:}), size(x));
end

function check_map(map, needed, name, sz)
    if needed
        assert(~isempty(map), 'apply_T2_LUT:MissingMap', ...
            'The LUT has a %s axis, so %smap is required.', name, name);
        assert(isequal(size(map), sz), 'apply_T2_LUT:MapSizeMismatch', ...
            '%smap size (%s) does not match apparentT2 size (%s).', name, mat2str(size(map)), mat2str(sz));
    else
        assert(isempty(map), 'apply_T2_LUT:UnusedMap', ...
            'The LUT has no %s axis, so %smap would be ignored.', name, name);
    end
end
