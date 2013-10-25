% == INITIALIZE ==
clear all;
close all;

% == INPUTS: DATA TO ANALYZE ==
% The directory with the groundtruth
sInputDir = '/home/thanuja/Dropbox/data/evaldata/';
sGroundTruth = 'labels2/';
% The directory with the experiments
sAlgoOutput = 'output/';
sExperiment = 'exp1/';

% == LOAD THE DATA == 
lDirGroundTruth = dir(fullfile(sInputDir,sGroundTruth));
lDirAlgoOutput = dir(fullfile(sInputDir,sAlgoOutput,sExperiment));

if(size(lDirGroundTruth) ~= size(lDirAlgoOutput))
    disp(['Error - You must provide directories with ', ...
          'the same corresponding number of file']);
      return;
else
    n=1;
    for i=1:length(lDirGroundTruth)
        if lDirGroundTruth(i).isdir == 0;
            disp(['Reading: groundtruth - ',lDirGroundTruth(i).name, ...
                  ', algooutput - ',lDirAlgoOutput(i).name]);
            
            cGroundTruth{n} = logical( ...
                                imread( ...
                                    fullfile( ...
                                        sInputDir, ...
                                        sGroundTruth, ... 
                                        lDirGroundTruth(i).name ...
                                            ) ...
                                       ) ...
                                      );
            cAlgoOutput{n} = logical( ...
                                imread( ...
                                    fullfile( ...
                                        sInputDir, ...
                                        sAlgoOutput, ...
                                        sExperiment, ...
                                        lDirAlgoOutput(i).name ...
                                            ) ...
                                       ) ...
                                      );
             n = n + 1;
        end
    end
end
 
% == CREATE CLUSTER IMAGES ==
[mGroundTruthClusters,nCC1] = createAllClusters(cGroundTruth);
[mAlgoOutputClusters,nCC2] = createAllClusters(cAlgoOutput);

figure('name','Clusters=[Label,Id]'),
     nL = max([nCC1,nCC2]);
     mCM = hsv2rgb([rand(nL,1),.7+0.3*rand(nL,1),.7+.3*rand(nL,1)]);
     
     subplot(1,2,1), colormap(polcmap(nCC1,0.8)), imagesc(mGroundTruthClusters), title('Ground Truth')
     subplot(1,2,2), colormap(polcmap(nCC2,0.8)), imagesc(mAlgoOutputClusters), title('Algo Output');
    
% == CLUSTERING BY LABEL & ID MEASURES ==
[nRI,nGCE,nVI] = compareSegmentations(double(mGroundTruthClusters), ...
                                      double(mAlgoOutputClusters));
disp(['Rand Index Clusters=[Label,Id]: ',num2str(nRI)]);
disp(['Variation of Information Clusters=[Label,Id]: ',num2str(nVI)]);
disp(['Global Consistency Error Clusters=[Label,Id]: ',num2str(nGCE)]);

% == CREATE LABEL IMAGES ==
mGroundTruthLabels = createAllLabels(cGroundTruth);
mAlgoOutputLabels = createAllLabels(cAlgoOutput);

figure('name','Clusters=[Label]'),
    subplot(1,2,1), colormap('gray'), imagesc(mGroundTruthLabels), title('Ground Truth')
    subplot(1,2,2), colormap('gray'), imagesc(mAlgoOutputLabels), title('Algo Output');
    
% == CLUSTERING BY LABEL MEASURES ==
[nRI,nGCE,nVI] = compareSegmentations(double(mGroundTruthLabels), ...
                                      double(mAlgoOutputLabels));
disp(['Rand Index Clusters=[Label]: ',num2str(nRI)]);
disp(['Variation of Information Clusters=[Label]: ',num2str(nVI)]);
disp(['Global Consistency Error Clusters=[Label]: ',num2str(nGCE)]);

% == PIXELWISE MEASURES==
[nPixAccTot,mPixMeasuresPerLabel] = comparePixelwise(double(mGroundTruthLabels), ...
                                                     double(mAlgoOutputLabels));

disp(['Total pixel accuracy: ',num2str(nPixAccTot)]);
disp('Per label pixel accuracy [Prec/Rec/Spec/NPV]:');
mPixMeasuresPerLabel


