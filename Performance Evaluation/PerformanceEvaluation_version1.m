%% ===================== Pixel-level Confusion Matrix from GT + Pred HDF masks (8-class, custom IDs) =====================
clc; clear; close all;

%% --------- USER SETTINGS (EDIT ONLY THESE PATHS) ---------
gtDir   = 'C:\Users\sriniva3\OneDrive - Aalto University\Litreature\Nokia Project\Spectrum Sensing\Wifi_vs_NR_baseline\baseline Spectum sensing DL\BaselineDetecionMethodPipeline_version2.0\dataset_out_180\gt8\maskhdf';      % e.g., ...\gt8\maskhdf
predDir = 'C:\Users\sriniva3\OneDrive - Aalto University\Litreature\Nokia Project\Spectrum Sensing\Wifi_vs_NR_baseline\baseline Spectum sensing DL\BaselineDetecionMethodPipeline_version2.0\dataset_out_180\baseline5\pred_maskhdf';    % e.g., ...\baseline8\pred_maskhdf  (or your pred folder)

% Your classes + pixel IDs (as provided)
classNames   = ["Noise" "NR" "LTE" "WLAN" "RADAR" "NRLTE" "NRWLAN" "LTEWLAN"];
pixelLabelID = [0, 63, 127, 191, 6, 32, 96, 160];

% Options
ignoreNoise = false;     % set true to ignore Noise pixels in evaluation
allowTransposeFix = true; % auto-fix GT vs Pred transpose mismatch
filePatterns = {'*.hdf','*.h5'}; % tries both

% Output folder
outDir = fullfile(predDir, 'confusion_eval_8class');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% --------- Collect files ---------
gtFiles = [];
predFiles = [];
for p = 1:numel(filePatterns)
    gtFiles   = [gtFiles;   dir(fullfile(gtDir,   filePatterns{p}))]; 
    predFiles = [predFiles; dir(fullfile(predDir, filePatterns{p}))]; 
end

if isempty(gtFiles),   error('No GT mask files found in: %s', gtDir); end
if isempty(predFiles), error('No Pred mask files found in: %s', predDir); end

gtMap   = buildIdMap_local(gtFiles);
predMap = buildIdMap_local(predFiles);

commonIds = intersect(gtMap.ids, predMap.ids);
if isempty(commonIds)
    error('No matching IDs between GT and Pred folders. Ensure filenames share the same numeric index.');
end
commonIds = sort(commonIds);

fprintf('Matched %d files.\n', numel(commonIds));

%% --------- Confusion accumulation ---------
K = numel(classNames);
C = zeros(K, K, 'double');   % rows=GT, cols=Pred

for k = 1:numel(commonIds)
    id = commonIds(k);

    gtPath   = fullfile(gtDir,   gtMap.byId(id));
    predPath = fullfile(predDir, predMap.byId(id));

    GT   = readMaskHDF_local(gtPath);
    PRED = readMaskHDF_local(predPath);

    GT   = double(GT);
    PRED = double(PRED);

    % Fix transpose mismatch if needed
    if ~isequal(size(GT), size(PRED))
        if allowTransposeFix && isequal(size(GT), fliplr(size(PRED)))
            PRED = PRED.';  % transpose prediction
        else
            warning('Size mismatch for ID=%d (GT %s vs Pred %s). Skipping.', id, mat2str(size(GT)), mat2str(size(PRED)));
            continue;
        end
    end

    gtVec   = GT(:);
    predVec = PRED(:);

    % Optional: ignore Noise (pixelLabelID(1) = 0)
    if ignoreNoise
        keep = (gtVec ~= pixelLabelID(1));
        gtVec = gtVec(keep);
        predVec = predVec(keep);
    end

    % Keep only known labels on both sides
    keep = ismember(gtVec, pixelLabelID) & ismember(predVec, pixelLabelID);
    gtVec = gtVec(keep);
    predVec = predVec(keep);

    % Map pixel IDs -> class indices (1..K)
    [~, gtIdx]   = ismember(gtVec, pixelLabelID);
    [~, predIdx] = ismember(predVec, pixelLabelID);

    % Accumulate
    C = C + accumarray([gtIdx, predIdx], 1, [K, K], @sum, 0);
end

%% --------- Metrics ---------
diagC  = diag(C);
rowSum = sum(C,2);
colSum = sum(C,1).';

precision = diagC ./ max(colSum, eps);
recall    = diagC ./ max(rowSum, eps);
f1        = 2*(precision.*recall) ./ max(precision+recall, eps);
iou       = diagC ./ max(rowSum + colSum - diagC, eps);

C_rowNorm = C ./ max(rowSum, eps);  % normalize rows (GT-normalized)

metricsTbl = table(classNames(:), pixelLabelID(:), precision, recall, f1, iou, rowSum, ...
    'VariableNames', {'Class','PixelID','Precision','Recall','F1','IoU','Support'});

disp('Confusion matrix (rows=GT class, cols=Pred class):');
disp(C);
disp(metricsTbl);

%% --------- Save outputs ---------
save(fullfile(outDir, 'confusion_results_8class.mat'), ...
    'C','C_rowNorm','classNames','pixelLabelID','ignoreNoise','commonIds','metricsTbl');

writetable(metricsTbl, fullfile(outDir, 'per_class_metrics.csv'));

% Confusion chart (row-normalized)
% Confusion chart (row-normalized; no extra summary row/col)
fig1 = figure('Color','w');
cc = confusionchart(C, cellstr(classNames));

cc.Normalization = 'row-normalized';
cc.RowSummary    = 'off';
cc.ColumnSummary = 'off';

cc.Title = 'Confusion Matrix-Synthetic(BaseLine based Detection ))';

saveas(fig1, fullfile(outDir, 'confusionchart_rowNorm.png'));

% % Raw counts heatmap
% fig2 = figure('Color','w');
% imagesc(C); axis image; colorbar;
% set(gca,'XTick',1:K,'XTickLabel',classNames, 'YTick',1:K,'YTickLabel',classNames);
% xtickangle(45);
% xlabel('Pred'); ylabel('GT');
% title('Confusion Matrix (raw pixel counts)');
% saveas(fig2, fullfile(outDir, 'confusion_counts.png'));

fprintf('Saved results to: %s\n', outDir);

%% ===================== Local helper functions =====================

function S = buildIdMap_local(files)
    ids = nan(numel(files),1);
    names = strings(numel(files),1);
    for i = 1:numel(files)
        names(i) = string(files(i).name);
        ids(i) = extractId_local(files(i).name);
    end
    keep = ~isnan(ids);
    ids = ids(keep);
    names = names(keep);

    [ids, order] = sort(ids);
    names = names(order);

    S.ids = ids;
    S.byId = containers.Map('KeyType','double','ValueType','char');
    for i = 1:numel(ids)
        S.byId(ids(i)) = char(names(i));
    end
end

function id = extractId_local(fname)
    % Extract the LAST digit-run in the filename: mask_0037.hdf -> 37
    tok = regexp(fname, '(\d+)(?!.*\d)', 'tokens', 'once');
    if isempty(tok), id = NaN; return; end
    id = str2double(tok{1});
end

function M = readMaskHDF_local(fpath)

    %% 1) Try imread first (HDF raster)
    try
        % Some HDF files are indexed-color images: [X,map] = imread(...)
        try
            [X, map] = imread(fpath);
            if ~isempty(map)
                % X already contains index values (often uint8). We want those indices/pixel values.
                M = X;
            else
                M = X;
            end
        catch
            % If the 2-output form fails, try 1-output
            M = imread(fpath);
        end

        % If RGB, take first channel (your masks should be single-channel IDs anyway)
        if ndims(M) == 3
            M = M(:,:,1);
        end

        M = double(squeeze(M));
        return;
    catch
        % fall through
    end

    %% 2) Try HDF5
    try
        info = h5info(fpath);
        candidates = {'/mask','/data','/labels','/label','/gt','/pred','/m','/dataset','/D'};
        for i = 1:numel(candidates)
            try
                M = h5read(fpath, candidates{i});
                M = double(squeeze(M));
                if ndims(M) > 2, M = M(:,:,1); end
                return;
            catch
            end
        end
        dsPath = firstDatasetPath_local(info);
        if ~isempty(dsPath)
            M = h5read(fpath, dsPath);
            M = double(squeeze(M));
            if ndims(M) > 2, M = M(:,:,1); end
            return;
        end
    catch
        % fall through
    end

    %% 3) Try HDF4 SDS/Vdata (if present)
    info4 = hdfinfo(fpath);

    if isfield(info4,'SDS') && ~isempty(info4.SDS)
        dsName = info4.SDS(1).Name;
        R = hdfread(fpath, dsName);
        M = unwrapHdfread_local(R);
        M = double(squeeze(M));
        if ndims(M) > 2, M = M(:,:,1); end
        return;
    end

    if isfield(info4,'Vdata') && ~isempty(info4.Vdata)
        vdName = info4.Vdata(1).Name;
        R = hdfread(fpath, vdName);
        M = unwrapHdfread_local(R);
        M = double(squeeze(M));
        if ndims(M) > 2, M = M(:,:,1); end
        return;
    end

    error('Could not read mask: %s (imread/h5read/hdfread all failed)', fpath);
end

function p = firstDatasetPath_local(info)
    p = '';
    if isfield(info,'Datasets') && ~isempty(info.Datasets)
        p = ['/' info.Datasets(1).Name];
        return;
    end
    if isfield(info,'Groups') && ~isempty(info.Groups)
        for g = 1:numel(info.Groups)
            gp = firstDatasetPath_local(info.Groups(g));
            if ~isempty(gp)
                p = [info.Groups(g).Name gp];
                return;
            end
        end
    end
end

function M = unwrapHdfread_local(R)
    M = R;
    if iscell(M)
        M = M{1};
        while iscell(M), M = M{1}; end
    end
    if isstruct(M)
        f = fieldnames(M);
        M = M.(f{1});
    end
end