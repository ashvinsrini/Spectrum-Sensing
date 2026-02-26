%% ===================== Pixel-level Confusion Matrix (ONLY NON-OVERLAP FILES + ONLY 5 BASE CLASSES) =====================
clc; clear; close all;

%% --------- USER SETTINGS (EDIT ONLY THESE PATHS) ---------
gtDir   = 'C:\Users\sriniva3\OneDrive - Aalto University\Litreature\Nokia Project\Spectrum Sensing\Wifi_vs_NR_baseline\baseline Spectum sensing DL\BaselineDetecionMethodPipeline_version2.0\dataset_out_180\gt8\maskhdf';
predDir = 'C:\Users\sriniva3\OneDrive - Aalto University\Litreature\Nokia Project\Spectrum Sensing\Wifi_vs_NR_baseline\baseline Spectum sensing DL\BaselineDetecionMethodPipeline_version2.0\dataset_out_180\baseline5\pred_maskhdf';

% ----- Keep ONLY these 5 classes in the confusion matrix -----
classNames   = ["Noise" "NR" "LTE" "WLAN" "RADAR"];
pixelLabelID = [0, 63, 127, 191, 6];   % Noise, NR, LTE, WLAN, RADAR

% Overlap IDs (used ONLY to exclude files whose GT contains any of these)
overlapIDs = [32, 96, 160];  % NRLTE, NRWLAN, LTEWLAN

% Options
ignoreNoise = false;        % if true, removes GT==Noise pixels from evaluation
allowTransposeFix = true;   % auto-fix GT vs Pred transpose mismatch
filePatterns = {'*.hdf','*.h5'}; % tries both

% Output folder
outDir = fullfile(predDir, 'confusion_eval_5class_NO_OVERLAP_FILES');
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
fprintf('Matched %d GT/Pred pairs.\n', numel(commonIds));

%% --------- Confusion accumulation (ONLY files with NO overlap in GT) ---------
K = numel(classNames);
C = zeros(K, K, 'double');   % rows=GT, cols=Pred

includedIds = [];
excludedIds = [];
nPredRemapped = 0;

for k = 1:numel(commonIds)
    id = commonIds(k);

    gtPath   = fullfile(gtDir,   gtMap.byId(id));
    predPath = fullfile(predDir, predMap.byId(id));

    % ---- Read GT first, exclude file if any overlap pixels exist ----
    GT = readMaskHDF_local(gtPath);
    GT = double(GT);

    if any(ismember(GT(:), overlapIDs))
        excludedIds(end+1,1) = id; 
        continue;
    end
    includedIds(end+1,1) = id; 

    % ---- Read PRED for included files ----
    PRED = readMaskHDF_local(predPath);
    PRED = double(PRED);

    % Fix transpose mismatch if needed
    if ~isequal(size(GT), size(PRED))
        if allowTransposeFix && isequal(size(GT), fliplr(size(PRED)))
            PRED = PRED.';  % transpose prediction
        else
            warning('Size mismatch for ID=%d (GT %s vs Pred %s). Skipping.', id, mat2str(size(GT)), mat2str(size(PRED)));
            includedIds(end) = []; % remove from included
            continue;
        end
    end

    gtVec   = GT(:);
    predVec = PRED(:);

    % Optional: ignore Noise pixels (GT-based)
    if ignoreNoise
        keep = (gtVec ~= pixelLabelID(1));
        gtVec = gtVec(keep);
        predVec = predVec(keep);
    end

    % Remap ANY prediction that is not in the 5-class set -> Noise (0)
    % This ensures a 5x5 confusion matrix without dropping those error pixels.
    badPred = ~ismember(predVec, pixelLabelID);
    nPredRemapped = nPredRemapped + nnz(badPred);
    predVec(badPred) = pixelLabelID(1); % Noise

    % Keep only GT pixels belonging to the 5-class set (should already be true if data is clean)
    keep = ismember(gtVec, pixelLabelID);
    gtVec = gtVec(keep);
    predVec = predVec(keep);

    % Map pixel IDs -> class indices (1..K)
    [~, gtIdx]   = ismember(gtVec,   pixelLabelID);
    [~, predIdx] = ismember(predVec, pixelLabelID);

    % Accumulate
    C = C + accumarray([gtIdx, predIdx], 1, [K, K], @sum, 0);
end

fprintf('Included (NO overlap in GT): %d files\n', numel(includedIds));
fprintf('Excluded (GT had overlap):   %d files\n', numel(excludedIds));
fprintf('Pred pixels remapped to Noise (pred not in 5-class set): %d\n', nPredRemapped);

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
save(fullfile(outDir, 'confusion_results_5class_NO_OVERLAP_FILES.mat'), ...
    'C','C_rowNorm','classNames','pixelLabelID','ignoreNoise','includedIds','excludedIds','metricsTbl','overlapIDs','nPredRemapped');

writetable(metricsTbl, fullfile(outDir, 'per_class_metrics.csv'));
writematrix(includedIds, fullfile(outDir, 'included_ids.txt'));
writematrix(excludedIds, fullfile(outDir, 'excluded_ids.txt'));

% Confusion chart (row-normalized; no extra summary row/col)
fig1 = figure('Color','w');
cc = confusionchart(C, cellstr(classNames));
cc.Normalization = 'row-normalized';
cc.RowSummary    = 'off';
cc.ColumnSummary = 'off';
cc.Title = 'Confusion Matrix-Synthetic(Baseline Detection)';
saveas(fig1, fullfile(outDir, 'confusionchart_rowNorm.png'));

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
    tok = regexp(fname, '(\d+)(?!.*\d)', 'tokens', 'once'); % last digit-run
    if isempty(tok), id = NaN; return; end
    id = str2double(tok{1});
end

function M = readMaskHDF_local(fpath)
% Robust: tries imread (incl. HDF raster via Ref), then HDF5, then HDF4 SDS/Vdata.

    % 1) Try imread directly
    try
        [X,map] = imread(fpath);
        if ~isempty(map), M = X; else, M = X; end
        if ndims(M)==3, M = M(:,:,1); end
        M = double(squeeze(M));
        return;
    catch
        % try HDF raster Ref
        try
            info = imfinfo(fpath);
            if isfield(info(1),'Ref')
                [X,map] = imread(fpath, 'Ref', info(1).Ref);
                if ~isempty(map), M = X; else, M = X; end
                if ndims(M)==3, M = M(:,:,1); end
                M = double(squeeze(M));
                return;
            end
        catch
        end
    end

    % 2) Try HDF5
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
    end

    % 3) Try HDF4 SDS/Vdata
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

    error('Could not read mask: %s (imread/imfinfo-Ref/h5read/hdfread all failed)', fpath);
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