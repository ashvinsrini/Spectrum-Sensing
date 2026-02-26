%% ===================== Stripe-level ROC for Baseline (from baseline_XXXX.mat) =====================
clc; clear; close all;
rng(1);

%% -------- USER: EDIT PATH --------
scoreDir = 'C:\Users\sriniva3\OneDrive - Aalto University\Litreature\Nokia Project\Spectrum Sensing\Wifi_vs_NR_baseline\baseline Spectum sensing DL\BaselineDetecionMethodPipeline_version2.0\dataset_out_180\baseline5\mat';  % where baseline_0001.mat etc live
outDir = fullfile(scoreDir, 'roc_eval_stripe');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% -------- Baseline score columns are for THESE 4 classes --------
% (Based on your baseline_0001.mat: stripeScores is 4x4 => 4 tech columns)
scoreClassOrder = ["NR","LTE","WLAN","RADAR"];   % assumed column order in stripeScores
K4 = numel(scoreClassOrder);

%% -------- Options --------
skipOverlapGT = true;   % ignore stripes whose GT is NRLTE/NRWLAN/LTEWLAN (baseline has no overlap score columns)
minValidScore = -inf;   % can set e.g. 0 if you know scores are nonnegative

%% -------- Load all baseline_*.mat --------
matFiles = dir(fullfile(scoreDir, 'baseline_*.mat'));
if isempty(matFiles), error('No baseline_*.mat found in %s', scoreDir); end

% per-class pooled labels/scores
allY = cell(K4,1);
allS = cell(K4,1);

for i = 1:numel(matFiles)
    P = fullfile(scoreDir, matFiles(i).name);
    B = load(P);

    if ~isfield(B,'stripeScores') || ~isfield(B,'labelsTrue')
        warning('Skipping %s (missing stripeScores/labelsTrue).', matFiles(i).name);
        continue;
    end

    scores = double(B.stripeScores);
    ytrue  = string(B.labelsTrue);
    ytrue  = ytrue(:);
    
    % --- FIX: enforce scores as Nstripes x K4 to match ytrue length ---
    scores = fixStripeScores_local(scores, numel(ytrue), K4);
    if isempty(scores)
        warning('Skipping %s (cannot reconcile stripeScores size with labelsTrue).', matFiles(i).name);
        continue;
    end
    
    % --- Now safe to build valid mask (same length as size(scores,1)) ---
    valid = ~ismissing(ytrue) & (strlength(ytrue) > 0);
    valid = valid & all(isfinite(scores),2) & any(scores > minValidScore, 2);
    
    scores = scores(valid,:);
    ytrue  = ytrue(valid);

    if isempty(ytrue), continue; end

    if skipOverlapGT
        isOverlap = (ytrue == "NRLTE") | (ytrue == "NRWLAN") | (ytrue == "LTEWLAN");
        scores = scores(~isOverlap,:);
        ytrue  = ytrue(~isOverlap);
    end

    % Build one-vs-rest pools
    for c = 1:K4
        allY{c} = [allY{c}; double(ytrue == scoreClassOrder(c))]; 
        allS{c} = [allS{c}; scores(:,c)]; 
    end
end

%% -------- Compute ROC + plot --------
auc = nan(K4,1);
roc = cell(K4,1);

fig = figure('Color','w'); hold on; grid on;
xlabel('False Positive Rate'); ylabel('True Positive Rate');
title('Baseline ROC (stripe-level, one-vs-rest)');

leg = strings(0);

for c = 1:K4
    y = allY{c};
    s = allS{c};

    if isempty(y) || numel(unique(y)) < 2
        warning('Class %s has insufficient positives/negatives. Skipping.', scoreClassOrder(c));
        continue;
    end

    % Prefer perfcurve if available; else use local ROC
    if exist('perfcurve','file') == 2
        [fpr,tpr,~,auc_c] = perfcurve(y, s, 1);
    else
        [fpr,tpr,auc_c] = rocCurve_local(y, s);
    end

    % If AUC < 0.5, try flipped sign (sometimes scores are "lower=better")
    if auc_c < 0.5
        if exist('perfcurve','file') == 2
            [fpr2,tpr2,~,auc2] = perfcurve(y, -s, 1);
        else
            [fpr2,tpr2,auc2] = rocCurve_local(y, -s);
        end
        if auc2 > auc_c
            fpr = fpr2; tpr = tpr2; auc_c = auc2;
            fprintf('Auto-flipped score polarity for %s (AUC improved).\n', scoreClassOrder(c));
        end
    end

    auc(c) = auc_c;
    roc{c} = struct('class',scoreClassOrder(c),'fpr',fpr,'tpr',tpr,'auc',auc_c);

    plot(fpr, tpr, 'LineWidth', 1.7);
    leg(end+1) = sprintf('%s (AUC=%.3f)', scoreClassOrder(c), auc_c); 
end

plot([0 1],[0 1],'k--'); % chance line
legend(leg, 'Location','southeast');

saveas(fig, fullfile(outDir, 'roc_stripe_level.png'));

resultsTbl = table(scoreClassOrder(:), auc, 'VariableNames', {'Class','AUC'});
disp(resultsTbl);
writetable(resultsTbl, fullfile(outDir, 'auc_stripe_level.csv'));
save(fullfile(outDir, 'roc_stripe_level.mat'), 'roc', 'resultsTbl', 'scoreClassOrder', 'skipOverlapGT');

fprintf('Saved ROC outputs to: %s\n', outDir);

%% -------- Local ROC (toolbox-free) --------
function [fpr,tpr,auc] = rocCurve_local(y, s)
    y = double(y(:));
    s = double(s(:));

    [~,ord] = sort(s,'descend');
    y = y(ord);

    P = sum(y==1); N = sum(y==0);
    if P==0 || N==0
        fpr=[0;1]; tpr=[0;1]; auc=NaN; return;
    end

    tp = cumsum(y==1);
    fp = cumsum(y==0);

    tpr = [0; tp./P];
    fpr = [0; fp./N];

    auc = trapz(fpr, tpr);
end


function Sfix = fixStripeScores_local(S, Nstr, K4)
% Make stripeScores into Nstripes x K4 to match labelsTrue length.
% Handles common layouts:
%   Nstr x K4   
%   K4 x Nstr  (transpose)
%   vector with Nstr*K4 elements (reshape)
% Returns [] if it cannot be reconciled.

    Sfix = [];

    if isempty(S) || ~isnumeric(S), return; end

    % squeeze just in case
    S = squeeze(S);

    if ismatrix(S)
        [r,c] = size(S);

        % already Nstr x K4
        if r == Nstr && c == K4
            Sfix = S; return;
        end

        % transposed: K4 x Nstr
        if r == K4 && c == Nstr
            Sfix = S.'; return;
        end

        % vector case: reshape
        if isvector(S) && numel(S) == Nstr*K4
            Sfix = reshape(S, [Nstr, K4]); return;
        end

        % Sometimes extra empty stripe row exists in one of them: allow trimming
        if c == K4 && r > Nstr
            Sfix = S(1:Nstr, :); return;
        end
        if c == K4 && r < Nstr
            % cannot invent missing rows reliably
            return;
        end
    end
end