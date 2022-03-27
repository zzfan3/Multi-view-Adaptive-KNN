function [kval]=FindMissingkvalue20151021(attrNode5,i,kval)
%% All samples of the left branch were used as the predictive sample of the right branch with no nodes
global attrNode;
kval=kval;
if ~isempty(attrNode(i).Lleaflabel)  % has left label
    kval=[kval;attrNode(i).kvalue];
    return;
end
if ~isempty(attrNode(i).Rleaflabel)  % has right label
    kval=[kval;attrNode(i).kvalue];
    return;
end
kval1=[];
kval2=[];
if ~isempty(attrNode(i).leftchildNode)&isempty(attrNode(i).Lleaflabel) 
    kval1=FindMissingkvalue20151021(attrNode,attrNode(i).leftchildNode,kval);

end
if ~isempty(attrNode(i).rightchildNode)&isempty(attrNode(i).Rleaflabel) 
    kval2=FindMissingkvalue20151021(attrNode,attrNode(i).rightchildNode,kval);

end
kval=[kval1;kval2];