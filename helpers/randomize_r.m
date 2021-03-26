function [pvalue, obs_stat, rand_stat, pvalue_corr] = randomize_r(rho, varargin)
% Sign test on accuracies.
% Inputs: correlation (length is number of observations/subjects).
% Name-value optional inputs: 'num_iterations' (default 5000) - number of randomizations;
%
% Outputs: p-value (one-tailed: number of randomizations exceeding observed statistic).
%          obs_stat (mean across observations)
%          rand_stat contains randomized statistic (length is num_iterations);
%
% I received this code from DC Dima (diana.c.dima@gmail.com)
p = inputParser;
addParameter(p, 'num_iterations',5000);
parse(p, varargin{:});

num_iterations = p.Results.num_iterations;
rand_sign = sign(randn(num_iterations,size(rho,1)));

if ismatrix(rho) %fewer than 3D, no memory issues
    if size(rho,1) == 1
        rho = rho(:);
        rand_rho = repmat(rho, 1, num_iterations,1).*rand_sign;
    elseif size(rho,2) == 1
        rand_rho = repmat(rho, 1, num_iterations)'.*rand_sign; %iterations x subjects
    elseif ~isvector(rho)
        %assume models are columns, subjects are rows
        num_models = size(rho,2);
        rand_rho = nan(num_iterations, size(rho,2),size(rho,1),'single'); %iterations x models x subjects
        for m = 1:num_models
            rand_rho(:,m,:) = repmat(rho(:,m), 1, num_iterations)'.*rand_sign;
        end
    end
    
    rand_stat = mean(rand_rho,ndims(rand_rho)); %average over subjects - the last dimension
    obs_stat = mean(rho,1); %average over subjects - the first dimension
    if isvector(rand_stat)
        pvalue = (length(find(rand_stat>=nanmean(rho)))+1)/(num_iterations+1);
    else
        pvalue = nan(1,num_models);
        pvalue_corr = nan(1,num_models);
        randmax = max(rand_stat,[],2);
        for m = 1:num_models
            pvalue(m) = (length(find(rand_stat(:,m)>=nanmean(rho(:,m))))+1)/(num_iterations+1);
            pvalue_corr(m) = (length(find(randmax>=nanmean(rho(:,m))))+1)/(num_iterations+1);
        end
        
    end
    
else
    
    if ndims(rho)==3 %sub x time x models
        
        sz1 = size(rho,2); sz2 = size(rho,3);
        rand_stat = nan(num_iterations,sz1,sz2);
        obs_stat = squeeze(nanmean(rho,1));
        for i = 1:num_iterations
            rsgn = repmat(squeeze(rand_sign(i,:))', [ 1 sz1 sz2]);
            rrho = rho.*rsgn;
            rand_stat(i,:,:) = nanmean(rrho,1);
        end
        
        pvalue = nan(size(obs_stat));
        for i = 1:size(obs_stat,1)
            for ii = 1:size(obs_stat,2)
                pvalue(i,ii) = (length(find(rand_stat(:,i,ii)>=obs_stat(i,ii)))+1)/(num_iterations+1);
            end
        end
        
        
    elseif ndims(rho)==4 %sub x space x time x models
        
        sz1 = size(rho,2); sz2 = size(rho,3); sz3 = size(rho,4);
        rand_stat = nan(num_iterations,sz1,sz2,sz3);
        obs_stat = squeeze(nanmean(rho,1));
        for i = 1:num_iterations
            rsgn = repmat(squeeze(rand_sign(i,:))', [ 1 sz1 sz2 sz3]);
            rrho = rho.*rsgn;
            rand_stat(i,:,:,:) = nanmean(rrho,1);
        end
        
        pvalue = nan(size(obs_stat));
        for i = 1:size(obs_stat,1)
            for ii = 1:size(obs_stat,2)
                for iii = 1:size(obs_stat,3)
                    pvalue(i,ii,iii) = (length(find(rand_stat(:,i,ii,iii)>=obs_stat(i,ii,iii)))+1)/(num_iterations+1);
                end
            end
        end

    end
        
end
