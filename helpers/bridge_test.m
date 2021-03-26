
function [b]=bridge_test(X1,X2,y, Lambda1, Lambda2)
% created by Haemy Lee Masson July/2020

Ytrain=y;
X1=normalize(X1);
X2=normalize(X2);
Ytrain=normalize(Ytrain);

[n1,f1] = size(X1); % n: observations, f: features
[n2,f2] = size(X2);

[n, v] = size(Ytrain); % n: observations, v: voxels
if n1~=n 
    error(message('stats:ridge:InputSizeMismatch')); 
end 
if n2~=n 
    error(message('stats:ridge:InputSizeMismatch')); 
end 

b=zeros(size(X1,2)+size(X2,2),size(y,2));
for v=1:size(y,2)
    lambda_one=Lambda1(v);
    lambda_two=Lambda2(v);
    bands = [ones(1,size(X1,2))*lambda_one ones(1,size(X2,2))*lambda_two];
    C = diag(bands);
    %Tikhonov banded ridge
    Xjoint = [X1 X2];
    LH = inv((Xjoint' * Xjoint) + (C' * C));
    XTY = Xjoint' * Ytrain(:,v);
    b(:,v) = LH * XTY;
end
end






