%% to obtain predictive class label
%% input: labelNode denotes class label belongs to which node; m denotes test samples numbers, Xte denotes test samples
%% output: predictLabel denotes predictive class label
function [predictLabel]=GetKTSlabel(labelNode,m,n,Xte)
global attrNode;
for i=1:m
    sample=attrNode(labelNode(i)).Sample(:,1:n);
    IND=KnnSearch(sample,Xte(i,:),1);  
    predictLabel(i)=attrNode(labelNode(i)).nearSampleLabel(IND);
end
predictLabel=predictLabel';