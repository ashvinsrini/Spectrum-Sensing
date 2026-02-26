imageSize = [256 256];                                 % pixels
sampleRate = 61.44e6;                                  % Hz
numSubFrames = 40;                                     % corresponds to 40 ms
frameDuration = numSubFrames*1e-3;                     % seconds
Dir_pth = '/home/sriniva3/spectrumsensing/RADAR/inputfiles_2026_onetechAndtwotechMix_per_spectrogram/';
trainDirRoot = Dir_pth;
trainingDataSource = "Downloaded data";
trainNow = false;
useCapturedData = false;
if trainingDataSource == "Generated data"
  numFramesPerStandard = 900;
  saveChannelInfo = false;
  helperSpecSenseTrainingData(numFramesPerStandard,classNames,imageSize, ...
      trainDirRoot,numSubFrames,sampleRate,saveChannelInfo);
end



baseNetwork = 'custom_LTE_NR_WLAN_RADAR_NoOverlap_2026';
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


classNames = ["Noise"   "NR"    "LTE"  "WLAN" "RADAR"];
pixelLabelID = [0, 63, 127, 191, 6];
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



%if ~strcmp(baseNetwork,"custom_LTE_NR_WLAN_RADAR_NoOverlap_2026")
    %layers = deeplabv3plus([256 256],numel(classNames),baseNetwork);
%end


%if strcmp(baseNetwork,"custom_LTE_NR_WLAN_RADAR_NoOverlap_2026")
    %layers = helperSpecSenseCustomNetwork(imageSize,numClasses);
%end


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
  MaxEpochs = 10, ...
  LearnRateSchedule = "piecewise",...
  InitialLearnRate = 0.02,...
  LearnRateDropPeriod = 10,...
  LearnRateDropFactor = 0.1,...
  ValidationData = cdsVal,...
  ValidationPatience = 5,...
  Shuffle="every-epoch",...
  OutputNetwork = "best-validation-loss",...
  Plots = 'none');





if trainNow
    [net,trainInfo] = trainnet(cdsTrain,layers, ...
        @(ypred,ytrue) lossFunction(ypred,ytrue,classWeights),opts); %#ok
    %save(sprintf('myNet_%s_%s',baseNetwork, ...
        %datetime('now',format='yyyy_MM_dd_HH_mm')), 'net')
        
        
        outdir = fullfile(pwd, 'models');

        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end 

        save(fullfile(outdir, 'model.mat'), 'net');
        fprintf('saving model\n');
        

else
    net = loadNetworkFromMATFile(baseNetwork); 
end


fprintf('training completed\n');



%%%%%%% test neural network %%%%%%%%%%%%%%%%
    testDir = '/home/sriniva3/spectrumsensing/RADAR/inputfiles_2026_onetechAndtwotechMix_per_spectrogram/';
    
    dataDir = testDir;
    imdsLTENRWLAN = imageDatastore(dataDir,FileExtensions=".png");
    pxdsResultsLTENRWLAN = semanticseg(imdsLTENRWLAN,net,MinibatchSize=mbs,WriteLocation=tempdir, ...
                            Classes=classNames); 



    pxdsTruthLTENRWLAN = pixelLabelDatastore(dataDir,classNames,pixelLabelID,...
                            FileExtensions=".hdf");
    metrics = evaluateSemanticSegmentation(pxdsResultsLTENRWLAN,pxdsTruthLTENRWLAN);



    %%%%% print normalized confusion matrix %%%%%%%%%%%%
    cm = confusionchart(metrics.ConfusionMatrix.Variables, ...
                        classNames, Normalization='row-normalized');
    cm.Title = 'Confusion Matrix - Deep Learning';



%%%%%%%%%%% plot received spectrogram %%%%%%%%%%%%%

  numSignals = length(imdsLTENRWLAN.Files);
  idx = 2
  rcvdSpectrogram1 = readimage(imdsLTENRWLAN,idx);
  trueLabels1 = readimage(pxdsTruthLTENRWLAN,idx);
  predictedLabels1 = readimage(pxdsResultsLTENRWLAN,idx);
  idx = length(imdsLTENRWLAN.Files);
  rcvdSpectrogram2 = readimage(imdsLTENRWLAN,idx);
  trueLabels2 = readimage(pxdsTruthLTENRWLAN,idx);
  predictedLabels2 = readimage(pxdsResultsLTENRWLAN,idx);

    figure(2);
    helperSpecSenseDisplayResults(rcvdSpectrogram1,trueLabels1,predictedLabels1, ...
      classNames,80e6,0,frameDuration)