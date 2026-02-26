% --- Select & reset the GPU MATLAB should use (respects CUDA_VISIBLE_DEVICES) ---
if canUseGPU
    g = gpuDevice([]);  % pick first visible GPU
    reset(g);
    disp(g);
else
    warning("GPU not available; training will fall back to CPU.");
end



file_pth = '/home/sriniva3/Documents/MATLAB/Examples/R2024b/deeplearning_shared/SpectrumSensingWithDeepLearning5GLTEExample';
imageSize = [256 256];                                 % pixels
sampleRate = 61.44e6;                                  % Hz
numSubFrames = 40;                                     % corresponds to 40 ms
frameDuration = numSubFrames*1e-3;                     % seconds
Dir_pth = '/home/sriniva3/spectrumsensing/RADAR/inputfiles_2026';
trainDirRoot = Dir_pth;
classNames = ["Noise"   "NR"    "LTE"  "WLAN" "RADAR"];
pixelLabelID = [0, 63, 127, 191, 6];
trainingDataSource = "Downloaded data";
trainNow = true;
useCapturedData = false;
if trainingDataSource == "Generated data"
  numFramesPerStandard = 900;
  saveChannelInfo = false;
  helperSpecSenseTrainingData(numFramesPerStandard,classNames,imageSize, ...
      trainDirRoot,numSubFrames,sampleRate,saveChannelInfo);
end



baseNetwork = 'custom';
trainDir = trainDirRoot; 
imageSize = [256 256];


%helperSpecSenseDownloadData(trainingDataSource,trainNow,useCapturedData, ...
    %baseNetwork,imageSize)

%%%%%%%%% load training data ############
folders = trainDir;
if useCapturedData
  folders = [folders,fullfile(trainDir,"captured")];
end

imds = imageDatastore(folders,FileExtensions=".png");

numClasses = length(classNames);

pxdsTruthLTENRWLAN = pixelLabelDatastore(folders,classNames,pixelLabelID,...
                                  FileExtensions=".hdf");



%%% Analyze data sets stats %%%%%%%%%%%%%%%

tbl = countEachLabel(pxdsTruthLTENRWLAN);
frequency = tbl.PixelCount/sum(tbl.PixelCount);
figure
bar(1:numel(classNames),frequency)
grid on
xticks(1:numel(classNames)) 
xticklabels(tbl.Name)
xtickangle(45)
ylabel("Frequency")
saveas(gcf, sprintf('plot_task%d.png', getenv('SLURM_ARRAY_TASK_ID')));
fprintf('Code completed until generation of the pixel historgam plot\n');

[imdsTrain,pxdsTrain,imdsVal,pxdsVal,imdsTest,pxdsTest] = ...
  helperSpecSensePartitionData(imds,pxdsTruthLTENRWLAN,[80 10 10]);
cdsTrain = combine(imdsTrain,pxdsTrain);
cdsVal = combine(imdsVal,pxdsVal);
cdsTest = combine(imdsTest,pxdsTest);



if ~strcmp(baseNetwork,"custom")
    %layers = deeplabv3plus([256 256],numel(classNames),baseNetwork);
    layers = helperSpecSenseCustomResNet18(imageSize,numClasses);
end


if strcmp(baseNetwork,"custom")
    layers = helperSpecSenseCustomNetwork(imageSize,numClasses);
end


imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
imageFreq(isnan(imageFreq)) = [];
classWeights = median(imageFreq) ./ imageFreq;
classWeights = classWeights/(sum(classWeights)+eps(class(classWeights)));
if length(classWeights) < numClasses
    classWeights = [classWeights; zeros(numClasses - length(classWeights),1)];
end


mbs = 40;
opts = trainingOptions("sgdm",...
  MiniBatchSize = mbs,...
  MaxEpochs = 2000, ...
  LearnRateSchedule = "piecewise",...
  InitialLearnRate = 0.001,...
  LearnRateDropPeriod = 1500,...
  LearnRateDropFactor = 0.1,...
  ValidationData = cdsVal,...
  ValidationPatience = 500,...
  Shuffle="every-epoch",...
  OutputNetwork = "best-validation-loss",...
  Plots = 'none',...
   ExecutionEnvironment = "gpu");





if trainNow
    [net,trainInfo] = trainnet(cdsTrain,layers, ...
        @(ypred,ytrue) lossFunction(ypred,ytrue,classWeights),opts); %#ok
    %save(sprintf('myNet_%s_%s',baseNetwork, ...
        %datetime('now',format='yyyy_MM_dd_HH_mm')), 'net')


        outdir = fullfile(pwd, 'models');

        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end 

        %save(fullfile(outdir, 'model_LTE_NR_WLAN_RADAR_new.mat'), 'net', 'trainInfo', '-v7.3');
        save(fullfile(outdir,'model_LTE_NR_WLAN_RADAR_new.mat'),'net','trainInfo');
        fprintf('saving model\n');


else
    net = loadNetworkFromMATFile(baseNetwork); 
end


fprintf('training completed\n');
