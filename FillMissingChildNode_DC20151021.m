function [attrNodeDC]=FillMissingChildNode_DC20151021(i,attrN)
global attrNodeDC;
global icount;
attrNo=attrN;
if isempty(attrNodeDC(i).rightchildNode)&isempty(attrNodeDC(i).Rleaflabel)&~isempty(attrNodeDC(i).leftchildNode)
    attrNodeDC(i).rightchildNode=icount+1;
    attrNodeDC(icount+1).id=icount+1;
    attrNodeDC(icount+1).Rleaflabel=10;  
    icount=icount+1;
end
if isempty(attrNodeDC(i).leftchildNode)&isempty(attrNodeDC(i).Lleaflabel)&~isempty(attrNodeDC(i).rightchildNode)
    attrNodeDC(i).leftchildNode=icount+1;
    attrNodeDC(icount+1).id=icount+1;
    attrNodeDC(icount+1).Lleaflabel=10;  
    icount=icount+1;
end
if ~isempty(attrNodeDC(i).Lleaflabel)|~isempty(attrNodeDC(i).Rleaflabel)
    return;
end
if ~isempty(attrNodeDC(i).leftchildNode)
    FillMissingChildNode_DC20151021(attrNodeDC(i).leftchildNode,attrNo);
end
if ~isempty(attrNodeDC(i).rightchildNode)
    FillMissingChildNode_DC20151021(attrNodeDC(i).rightchildNode,attrNo);
end 
    