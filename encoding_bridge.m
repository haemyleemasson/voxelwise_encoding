% voxel-wise encoding model with banded ridge regression
% created by haemy lee masson July/2020
% dependency: cosmo mvpa toolbox
clearvars; clc; close all;
addpath(fullfile(pwd,'helpers'));
cd ../../;
basepath = pwd;
fmripath = fullfile(basepath, 'fmriData/');

model_name={'full'}
model=model_name{m}
disp(['model: ', model])
groupdir = fullfile(fmripath, ['results/group/' model '/']);
if ~exist(groupdir, 'dir')
    mkdir(groupdir);
end
featuredir = fullfile(basepath, 'analysis/features/'); % where are the features?
name_feature={'hue', 'saturation', 'pixel', 'valence', 'arousal'}; % name your features
load('feature1.mat') %load your feature space ; a feature space that use the same lambda
load('feature2.mat') %load another feature space (to use banded ridge, you may want to provide two or more feature spaces).
nfeature=size(name_feature,2);
X1=normalize(feature1); 
X2=normalize(feature2);
%% subject level analysis, voxel-wise encoding with ridge regression.
alphas = logspace(1,4,10);
ratios = logspace(-2,2,15);
for subj=1
    tic
    disp(['subject: ' num2str(subj)])
    clear data_clean;
    savedir = fullfile(fmripath, ['results/sub' num2str(subj) '/' model '/']);
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
    [data_clean] = clean_image(filename, subj, mask, fmripath); %load image, mask it, and clean the data (Nan and Inf)
    numvoxel=size(data_clean.samples,2);  %how many voxels?
    y=data_clean.samples; %BOLD signal
    tr=size(y,1); %how many TRs?
    niteration=10;
    nfold=5; %how many folds?
    r_iterated=zeros(numvoxel,niteration); %define the variable, r
    weights_interated=zeros(nfeature, numvoxel, niteration); %define the variable, weights
    b_interated=zeros(size([X1 X2],2), numvoxel, niteration); %define the variable, weights
    c = cvpartition(tr,'KFold',niteration);
    for iter=1:c.NumTestSets %iteration k fold
        idxTrain = training(c,iter); %indices of training set
        idxTest = ~idxTrain; % indices of test set
        Lmat1 = zeros(numvoxel,length(alphas)*length(ratios), nfold,'single');
        Lmat2 = zeros(numvoxel,length(alphas)*length(ratios), nfold,'single');
        %define the variable Lmat (number of voxels * number of lambdas * number of folds)
        X1_training=X1(idxTrain,:);
        X2_training=X2(idxTrain,:);
        y_training=y(idxTrain,:);
        innerc = cvpartition(c.TrainSize(iter),'KFold',nfold); % Kfold cross-validation for the hyper-parameter
        for fold=1:nfold % find the best lambda across kfolds
            disp(['itertation: ', num2str(iter), ' fold: ', num2str(fold)]);
            SSE1_all=zeros(numvoxel,length(ratios),length(alphas),'single');
            SSE2_all=zeros(numvoxel,length(ratios),length(alphas),'single');
            innerTrain = training(innerc,fold); %indices of training set
            innerTest = ~innerTrain; %indices of validation set for lambda
            [b_banded, lamb1, lamb2]=bridge(X1_training(innerTrain,:),X2_training(innerTrain,:), y_training(innerTrain,:),alphas, ratios, 3);
            disp('training done in the inner loop');
            SSE1=zeros(numvoxel,length(ratios),'single');
            SSE2=zeros(numvoxel,length(ratios),'single');
            for i=1:size(b_banded,4)
                for k=1:size(b_banded,3)
                    innerX1=X1_training(innerTest,:);
                    innerX1(isnan(innerX1)) = 0;
                    innerX2=X2_training(innerTest,:);
                    innerX2(isnan(innerX2)) = 0;
                    y1_hat = innerX1 * b_banded(1:size(innerX1,2),k,i); %predicted bold signal for each lambda
                    y2_hat = innerX2 * b_banded(size(innerX1,2)+1:end,k,i); %predicted bold signal for each lambda
                    y_true = y_training(innerTest,:); % ground truth
                    error1 = (y_true-y1_hat); %error
                    error2 = (y_true-y2_hat); %error
                    SSE1(:,k) = sum(error1.^2); %SSE for each lambda
                    SSE2(:,k) = sum(error2.^2); %SSE for each lambda
                end
                SSE1_all(:,:,i) = SSE1;
                SSE2_all(:,:,i) = SSE2;
            end
            SSE1_all_reshape = reshape(SSE1_all,[size(SSE1_all,1),size(SSE1_all,2)*size(SSE1_all,3)]);
            SSE2_all_reshape = reshape(SSE2_all,[size(SSE2_all,1),size(SSE2_all,2)*size(SSE2_all,3)]);
            Lmat1(:,:,fold) = SSE1_all_reshape; %store SSE for each fold and each lambda
            Lmat2(:,:,fold) = SSE2_all_reshape;
            disp('testing done in the inner loop');
        end
        [~, Lambda1] = min(mean(Lmat1,3),[],2); % give me the best lambda with the smallest SSE
        [~, Lambda2] = min(mean(Lmat2,3),[],2);
        Lambda1 = lamb1(Lambda1); %lambda is selected for each voxel
        Lambda2 = lamb2(Lambda2);
        disp('selecting best hyperparameters done');
        [b]=bridge_test(X1(idxTrain,:),X2(idxTrain,:), y(idxTrain,:),Lambda1, Lambda2);
        disp('training for testing done');
        X=[X1 X2];
        outerX=X(idxTest,:);
        outerX(isnan(outerX)) = 0;
        y_true = y(idxTest,:); % ground truth
        [R, R_features] = prediction(b,outerX,y_true,name_feature,featuredir);
        r_iterated(:,iter)=R'; %store the R from the each iteration
        weights_interated(:, :, iter)=R_features; %store the beta weights from the each iteration
        b_interated(:, :, iter)=b;
        disp('testing done');
    end
    r_mean=mean(r_iterated,2)';
    weight_mean=mean(weights_interated,3); %R value for each feature space.
    b_mean=mean(b_interated,3);
    % save the results and make the nii file for each subject and weight
    [group_r{subj}, group_weights{subj}, group_weights_all{subj}] = bridge_results_nii(model, data_clean, subj, r_mean, b_mean, weight_mean, savedir, name_feature);
    save([groupdir 'group_r.mat'], 'group_r'); %r for the model
    save([groupdir 'group_weights.mat'], 'group_weights'); % r for each feature
    save([groupdir 'group_weights_all.mat'], 'group_weights_all'); % beta weights for all 
    disp(['duration: ', num2str(round(toc/60)), ' mins'])
end



