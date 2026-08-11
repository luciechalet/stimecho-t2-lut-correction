function [correctedT2, info] = apply_T2_LUT(apparentT2, LUT, args)
    arguments 
        apparentT2 {mustBeNumeric}
        LUT (1,1) struct
        args.mask {mustBeNumericOrLogical} = []
    end

    % Apply the T2 LUT to correct stimulated echoes in apparent T2 values.
    
    % Usage:
    %  [correctedT2, info] = apply_T2_LUT(apparentT2, LUT);
    %  [correctedT2, info] = apply_T2_LUT(apparentT2, LUT, 'mask', ROI_mask);

    % Inputs:
    %  apparentT2: Apparent T2 values array (must be in ms) 
    %  LUT: Look up table structure built with build_T2_LUT
    %  args.mask: Optional logical or numeric mask to apply the correction only to selected voxels (voxels outside the mask will be returned as NaN)

    % Outputs:
    %  correctedT2: Corrected T2 values array (in ms) same dimensions as apparentT2
    %  info: Structure containing information about the LUT and the correction process

    %% Sanity checks: LUT structure
    required_fields = {'TrueT2', 'ApparentT2', 'T2Units'};
    missing_fields = required_fields(~isfield(LUT, required_fields));
    assert(isempty(missing_fields), ...
        'LUT is missing field(s): %s. Expected a struct from jsondecode of a build_T2_LUT JSON.', ...
        strjoin(missing_fields, ', '));
    assert(strcmpi(LUT.T2Units, 'ms'), ...
        'LUT.T2Units is %s, but apply_T2_LUT expects LUT in ms.', LUT.T2Units);

    apparentAxis = LUT.ApparentT2(:);
    trueAxis     = LUT.TrueT2(:);
    assert(numel(apparentAxis) == numel(trueAxis), ...
        'LUT.ApparentT2 (%d) and LUT.TrueT2 (%d) differ in length.', ...
        numel(apparentAxis), numel(trueAxis));
    assert(all(diff(apparentAxis) > 0), ...
        'LUT.ApparentT2 is not strictly increasing, so the table cannot be inverted.');

    %% Apply the LUT to correct apparent T2 values
    x = double(apparentT2); % ensure double for interp1
    correctedT2 = interp1(apparentAxis, trueAxis, x, 'pchip', NaN); % NaN for out-of-bounds values, pchip for shape preserving cubic interpolation

    %% Apply mask if provided
    if ~isempty(args.mask)
        assert(isequal(size(args.mask), size(apparentT2)), ...
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
        'RangeLUT', [apparentAxis(1), apparentAxis(end)], ...
        'Units', 'ms');
end