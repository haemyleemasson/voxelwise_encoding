function [group_r,group_weights,group_weights_all] = bridge_results_nii(model, data_clean, subj,r_mean, b_mean, weight_mean,savedir,name_feature)
%save data to nii format
data_clean.samples=r_mean;
encoding_all=data_clean;
cosmo_map2fmri(encoding_all,[savedir model '_r_sub' num2str(subj) '.nii']);
group_r=encoding_all;
group_weights_all=b_mean;
for i=1:size(name_feature,2)
    data_clean.samples=weight_mean(i,:);
    encoding_weights=data_clean;
    cosmo_map2fmri(encoding_weights,[savedir name_feature{i} '_sub' num2str(subj) '.nii']);
    save_weights(i,:)=data_clean.samples;
end
data_clean.samples=save_weights;
group_weights=data_clean;
end



