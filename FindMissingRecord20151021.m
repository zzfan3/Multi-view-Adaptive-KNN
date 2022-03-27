function [record]=FindMissingRecord20151021(attrNode4,i,record)
%% All samples of the left branch were used as the predictive sample of the right branch with no nodes
global attrNode;
record=record;
if ~isempty(attrNode(i).Lleaflabel)  % has left label
    record=[record;attrNode(i).Record];
    return;
end
if ~isempty(attrNode(i).Rleaflabel)  % has right label
    record=[record;attrNode(i).Record];
    return;
end
record1=[];
record2=[];
if ~isempty(attrNode(i).leftchildNode)&isempty(attrNode(i).Lleaflabel) 
    record1=FindMissingRecord20151021(attrNode,attrNode(i).leftchildNode,record);

end
if ~isempty(attrNode(i).rightchildNode)&isempty(attrNode(i).Rleaflabel) 
    record2=FindMissingRecord20151021(attrNode,attrNode(i).rightchildNode,record);

end
record=[record1;record2];