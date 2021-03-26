
function [b_banded, lamb1, lamb2]=bridge(X1_training,X2_training, y_training,alphas, ratios, method)
% Modifed Matlab version of Tikhonov regression function (tikreg) in Python
% with a polar grid search using ratio and alpha
% method 1 -> standard Tikhonov regression; method 2-> standard form solution; 
% method 3 -> Singular value decomposition (the fastest solution!)
% created by Haemy Lee Masson July/2020
% Citation: Nunez-Elizalde AO, Huth AG, and Gallant, JL (2019). 
% Voxelwise encoding models with non-spherical multivariate normal priors. NeuroImage.

X1=normalize(X1_training); % z-score: mean=0, std=1
X2=normalize(X2_training);
Ytrain=normalize(y_training);

[n1,f1] = size(X1); % n: observations, f: features
[n2,f2] = size(X2);

[n, v] = size(Ytrain); % n: observations, v: voxels
if n1~=n 
    error(message('stats:ridge:InputSizeMismatch')); 
end 
if n2~=n 
    error(message('stats:ridge:InputSizeMismatch')); 
end 
    
angle = atan(ratios); % ratio to angle for lamda polar grid search
b_banded=zeros(size(X1,2)+size(X2,2),size(Ytrain,2),length(angle),length(alphas)); % size(f1+f2, voxels, angles, scaling)
lamb1=zeros(length(angle),length(alphas)); % angle and alpha (scaling) determine lamba1 and 2 in polar grid search
lamb2=zeros(length(angle),length(alphas)); 
for k=1:length(alphas)
    solution=zeros(size(X1,2)+size(X2,2),size(Ytrain,2),length(angle)); % beta
    lambda1=zeros(length(angle),1); 
    lambda2=zeros(length(angle),1);
    for i=1:length(angle)
        lambda_one = cos(angle(i))*alphas(k); %lambda1 
        lambda_two = sin(angle(i))*alphas(k); %lambda2 
        bands = [ones(1,size(X1,2))*lambda_one ones(1,size(X2,2))*lambda_two]; 
        C = diag(bands);
        Cinv = inv(C); % C is a diagonal matrix with lambda1 and 2. It is used to panalize beta weights. 
        if method == 1
            %Tikhonov banded ridge (computationally very expensive)
            Xjoint = [X1 X2];
            LH = inv((Xjoint' * Xjoint) + (C' * C));
            XTY = Xjoint' * Ytrain;
            solution(:,:,i) = LH * XTY;
        elseif method == 2
            %Standard scaled banded ridge  (computationally less expensive)
            A = [X1/(lambda_one/alphas(k)) X2/(lambda_two/alphas(k))];
            LH = inv((A' * A) + ((alphas(k)^2) * eye(size(A,2))));
            RH = A' * Ytrain;
            solution_standard = (LH * RH) * alphas(k);
            solution(:,:,i) = Cinv * solution_standard; 
        else
            % Banded ridge with SVD (fastest)
            A = [X1/(lambda_one/alphas(k)) X2/(lambda_two/alphas(k))];
            [U, S, V] = svd(A,'econ');
            UTY = U' * Ytrain;
            D = diag(S ./ (S^2 + alphas(k)^2));
            D = diag(D);
            solution_svd = (V * D * UTY) * alphas(k); 
            solution(:,:,i) = Cinv * solution_svd; 
        end
    lambda1(i)=lambda_one;
    lambda2(i)=lambda_two;
    end
    b_banded(:,:,:,k)=solution; % beta value for each feature, voxel, and lambda
    lamb1(:,k)=lambda1; % lambda1
    lamb2(:,k)=lambda2; % lambda2
end
end





