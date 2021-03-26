clearvars; clc; close all;
addpath(fullfile(pwd,'helpers'));
cd ../../;
basepath = pwd;
groupdir = fullfile(basepath, 'fmriData/results/group/');
feature={'full', 'full-social', 'full-perceptual'}; 

for i=1:length(feature)
    group_corr=load([groupdir feature{i} '/group_r.mat'])
    group=group_corr.group_r;
for numsub=1:size(group,2)
    group{1,numsub}.samples=(group{1,numsub}.samples).^2;
    mean_group_r(numsub)=sqrt(mean(group{1,numsub}.samples));
end
mean(mean_group_r)
[idxs,group_intersect_cell]=cosmo_mask_dim_intersect(group); % remove un-shared voxels across subj
nsubj=numel(group_intersect_cell);
for subject_i=1:nsubj %tedious job to do.. we need to put chunks and target info for ttest.
    stacked_group=group_intersect_cell{subject_i};
    stacked_group.samples=group_intersect_cell{subject_i}.samples;
    stacked_group.sa.chunks=ones(1,1)*subject_i; %number of subj
    stacked_group.sa.targets=ones(1,1); 
    group_intersect{subject_i}=stacked_group;
end
group_clean=cosmo_stack(group_intersect,1,'drop_nonunique');
group_data{1,i}=group_clean.samples;
end

social=(group_data{1}-group_data{2});
perceptual=(group_data{1}-group_data{3});
[p_social, obs_stat_social, rand_stat_social, pvalue_corr] = randomize_r(social);
[p_perceptual, obs_stat_perceptual, rand_stat_perceptual, pvalue_corr] = randomize_r(perceptual);
[h_social, crit_p, adj_ci_cvrg, adj_p_social]=fdr_bh(p_social,.05);
[h_perceptual, crit_p, adj_ci_cvrg, adj_p_perceptual]=fdr_bh(p_perceptual,.05);

group_social=mean(social);
group_perceptual=mean(perceptual);
group_social(find(h_social==0))=0;
group_perceptual(find(h_perceptual==0))=0;
group_clean.sa.chunks=ones(1,1)*1; %number of subj
group_clean.sa.targets=ones(1,1); 
group_clean.samples=group_social;
cosmo_map2fmri(group_clean,[groupdir 'group_permu_unique_social.nii']);
group_clean.samples=group_perceptual;
cosmo_map2fmri(group_clean,[groupdir 'group_permu_unique_perceptual.nii']);

common=(group_data{2}+group_data{3}-group_data{1});
[p_common, obs_stat_common, rand_stat_common, pvalue_corr] = randomize_r(common);
[h_common, crit_p, adj_ci_cvrg, adj_p_common]=fdr_bh(p_common,.05);
group_common=mean(common);
group_common(find(h_common==0))=0;
group_clean.samples=group_common;
cosmo_map2fmri(group_clean,[groupdir 'group_permu_common_perceptual_social.nii']);

