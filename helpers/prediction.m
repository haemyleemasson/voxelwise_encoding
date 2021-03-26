function [R, R_features] = prediction(b,outerX,y_true,name_feature,featuredir)
% correlation between estimated y and true y.
% made by haemy lee masson July/2020
numvoxel=size(b,2);
y_hat = outerX * b; %predicted bold signal
R=zeros(numvoxel,1,'single');
for v=1:numvoxel
     R(v) = corr(y_true(:,v), y_hat(:,v));
end
k=1;
R_features=zeros(size(name_feature,2),numvoxel);
weight=zeros(numvoxel,1,'single');
for i=1:size(name_feature,2)
    feature=[featuredir name_feature{i}, '.mat'];
    A = importdata(feature);
    if size(A,1)*size(A,2)>1976
        y_hat = outerX(:,k:k+size(A,2)-1) * b(k:k+size(A,2)-1,:);
        k=k+((size(A,1)*size(A,2))/1921);
    else
        y_hat = outerX(:,k) * b(k,:);
        k=k+((size(A,1)*size(A,2))/1976);
    end
    for v=1:numvoxel
       weight(v) = corr(y_true(:,v), y_hat(:,v));
    end
    R_features(i,:)=weight';
end
end

