function require_stimfit(stimfit_dir)
    % Put StimFit on the MATLAB path. 
    % Resolution order:
    %     1. explicit argument, if given
    %     2. already on the path      -> nothing to do
    %     3. stimfit_path.txt in the repo root (local override, gitignored)
    %     4. sibling folder ../StimFit
    %     5. STIMFIT_DIR env var, if someone happens to have one
    %     6. ask once, and remember the answer
    % Usage:
    %       require_stimfit(); % add StimFit to path if not already there
    % Inputs:
    %       stimfit_dir: optional path to StimFit folder (default: ask user if not

    here   = fileparts(fileparts(mfilename('fullpath')));   % private/ -> repo root
    cfg    = fullfile(here, 'stimfit_path.txt');
    parent = fileparts(here);
    asked  = false;

    if nargin > 0 && ~isempty(stimfit_dir)
        d = stimfit_dir;
    elseif exist('StimFit_optset','file') == 2
        return
    elseif exist(cfg,'file') == 2
        d = strtrim(fileread(cfg));
    elseif exist(fullfile(parent,'StimFit','Engine','StimFit_optset.m'),'file') == 2
        d = fullfile(parent,'StimFit');
    elseif ~isempty(getenv('STIMFIT_DIR'))
        d = getenv('STIMFIT_DIR');
    else
        d = uigetdir(parent, 'Locate the StimFit toolbox folder');
        if isequal(d,0)
            error('require_stimfit:Cancelled', ...
                'StimFit is required but was not located.');
        end
        asked = true;
    end

    if exist(d,'dir') ~= 7
        error('require_stimfit:BadPath','StimFit folder does not exist: %s', d);
    end
    addpath(genpath(d));
    if exist('StimFit_optset','file') ~= 2
        error('require_stimfit:NotStimFit', ...
            'StimFit_optset not found under %s - wrong folder?', d);
    end

    if asked                     
        fid = fopen(cfg,'w');
        if fid > 0
            fprintf(fid,'%s\n',d);
            fclose(fid);
            fprintf('Remembered StimFit location in %s\n', cfg);
        end
    end
end