clearvars; clc; close all;
addpath(fullfile(pwd,'helpers'));
cd ../../;
basepath = pwd;
groupdir = fullfile(basepath, 'fmriData/results/new_group/full/'); % where the final results are? change the folder name
name_feature={'hue','saturation','pixel','MotionEnergy', 'text','indoor','face',...
    'amplitude','pitchHZ','music', 'sociality', 'speaking',...
    'mentalization','valence', 'arousal', 'layer5', }; % a list of features. change here based on your need. 
load([groupdir 'group_weights.mat']) % this file contains beta values for all featuers. 
for numsub=1:size(group_weights,2)
    group_weights{1,numsub}.samples=atanh(group_weights{1,numsub}.samples);
end
[idxs,group_intersect_cell]=cosmo_mask_dim_intersect(group_weights); % remove un-shared voxels across subj
nsubj=numel(group_intersect_cell);
for subject_i=1:nsubj %tedious job to do.. we need to put chunks and target info
    stacked_group=group_intersect_cell{subject_i};
    stacked_group.samples=group_intersect_cell{subject_i}.samples;
    group_intersect{subject_i}=stacked_group;
    feature_r{subject_i}=group_intersect{1,subject_i}.samples;
end
group_clean=cosmo_stack(group_intersect,1,'drop_nonunique');
feature=reshape(cell2mat(feature_r), [size(group_intersect{1,1}.samples,1), size(group_intersect{1,1}.samples,2), size(feature_r,2)]) ;
feature=feature(1:size(name_feature,2),:,:);
mean_feature=median(feature,3);
imagesc(mean_feature)
yticks([1:size(name_feature,2)])
yticklabels(name_feature);
xlabel('voxels')
median(mean_feature,2)

[~, ind1] = max(mean(feature,3),[],1);
group_clean.samples=ind1;
cosmo_map2fmri(group_clean,[groupdir 'group_feature_prefer_map.nii']);
for i=1:size(ind1,2)
mean_feature(ind1(1,i),i)=-inf;
end
[max2, ind2] = max(mean_feature);
group_clean.samples=ind2;
cosmo_map2fmri(group_clean,[groupdir 'group_feature_prefer_second_map.nii']);

