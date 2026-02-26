if canUseGPU
    g = gpuDevice([]);  % pick first visible GPU
    reset(g);
    disp(g);
else
    warning("GPU not available; training will fall back to CPU.");
end



%file_pth = 'C:\Users\sriniva3\OneDrive - Aalto University\Documents\MATLAB\Examples\R2025b\deeplearning_shared\SpectrumSensingWithDeepLearning5GLTEExample';
file_pth = 'C:\Users\sriniva3\OneDrive - Aalto University\Litreature\Nokia Project\Spectrum Sensing\Wifi_vs_NR_baseline\baseline Spectum sensing DL\Data Creation\dataset_out_mix_withoverlap\inputfiles_2026';
imageSize = [256 256];                                 % pixels
sampleRate = 61.44e6;                                  % Hz
numSubFrames = 40;                                     % corresponds to 40 ms
frameDuration = numSubFrames*1e-3;                     % seconds
trainDirRoot = fullfile(file_pth);
classNames = ["Noise" "NR" "LTE" "WLAN" "RADAR" "NRLTE" "NRWLAN" "LTEWLAN"];
%classNames = ["Noise" "NR" "LTE" "WLAN" "RADAR" ];
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
trainDir = fullfile(trainDirRoot); 
imageSize = [256 256];


%helperSpecSenseDownloadData(trainingDataSource,trainNow,useCapturedData, ...
    %baseNetwork,imageSize)

%%%%%%%%% load training data ############baseNetwork
folders = trainDir;
if useCapturedData
  folders = [folders,fullfile(trainDir,"captured")];
end
imds = imageDatastore(folders,FileExtensions=".png");


numClasses = length(classNames);
%pixelLabelID = floor((0:numClasses-1)/(numClasses-1)*255);
%pixelLabelID = [pixelLabelID(1:end-1),6];
pixelLabelID = [0, 63, 127, 191, 6, 32, 96, 160];
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


%helperSpecSenseDownloadData(trainingDataSource,trainNow,useCapturedData, ...
    %baseNetwork,imageSize)