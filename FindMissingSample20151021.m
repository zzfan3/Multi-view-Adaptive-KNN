function [samp]=FindMissingSample20151021(AttrNode,i,samp)
%% All samples of the left branch were used as the predictive sample of the right branch with no nodes
samp=samp;
if ~isempty(AttrNode(i).Lleaflabel)  % has left label
    samp=[samp;AttrNode(i).Sample];
    return;
end
if ~isempty(AttrNode(i).Rleaflabel)  % has right label
    samp=[samp;AttrNode(i).Sample];
    return;
end
fsample1=[];
fsample2=[];
if ~isempty(AttrNode(i).leftchildNode)&isempty(attrNode(i).Lleaflabel) % 修改过
    fsample1=FindMissingSample20151021(AttrNode,AttrNode(i).leftchildNode,samp);

end
if ~isempty(AttrNode(i).rightchildNode)&isempty(attrNode(i).Rleaflabel) % 修改过
    fsample2=FindMissingSample20151021(AttrNode,AttrNode(i).rightchildNode,samp);

end
samp=[fsample1;fsample2];