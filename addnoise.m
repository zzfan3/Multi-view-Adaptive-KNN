function [X]=addnoise(X)

viewnum = length(X);

for i=1:viewnum
    randIdxtr=randperm(size(X{i},1));
    noisenum=fix(0.20*size(X{i},1));
    XN{i} = X{i}(randIdxtr(1:noisenum),:);
    noisy=randn(size(XN{i}));
    X{i}(randIdxtr(1:noisenum),:)=noisy+X{i}(randIdxtr(1:noisenum),:);
end