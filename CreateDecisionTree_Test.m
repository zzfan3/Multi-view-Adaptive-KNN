%% build decision tree, don't perserve leaf node samples (ID3)
%% attrNode denotes built tree
%% D denotes training samples + labels
%% attrlist is used for record weather the atrribute is use, initial values are [1 2 3 4 ....], denotes [attribute1, attribute2, attrbute3,...]
%% i denotes generated i-th node
%% iparenttag is used for record father node variable
%% j=1 denotes generate left child node, j=2 denotes generate right child node
function [attrNodeDC]=CreateDecisionTree_Test(D,attrlist,i,iparenttag,j)
global attrNodeDC;
global icount;
[r, c]=size(D);
X=D(:,1:c-1);
Y=D(:,c);
attrNodeDC(i).id=i;
icount=icount+1;
if length(unique(Y))==1      % All the same class
    if j==1
        attrNodeDC(i).Lleaflabel=Y(1);
        attrNodeDC(iparenttag).leftchildNode=attrNodeDC(i).id;
    end
    if j==2
        attrNodeDC(i).Rleaflabel=Y(1);
        attrNodeDC(iparenttag).rightchildNode=attrNodeDC(i).id;
    end
    return;
end
if isempty(attrlist)     
    uy=unique(Y);
    fy=histc(Y,uy);
    [mv mi]=max(fy);
    if j==1
        attrNodeDC(i).Lleaflabel=uy(mi);
        attrNodeDC(iparenttag).leftchildNode=attrNodeDC(i).id;
    end
    if j==2
        attrNodeDC(i).Rleaflabel=uy(mi);
        attrNodeDC(iparenttag).rightchildNode=attrNodeDC(i).id;
    end
    return;
end
[splitattr, splitpoint]=AttributeSelectMethod(D,attrlist);  %Through the attribute selection method to select the split attribute, split point
attrNodeDC(i).splitattr=splitattr;
attrNodeDC(i).splitpoint=splitpoint;
if j==1                 %  The current node is the left child of his father's node
    attrNodeDC(iparenttag).leftchildNode=attrNodeDC(i).id;
end
if j==2
    attrNodeDC(iparenttag).rightchildNode=attrNodeDC(i).id;
end
iparenttag=i;
t=1;
for al=1:length(attrlist)  % 
    if attrlist(al)~=splitattr
        attrlistnew(t)=attrlist(al);
        t=t+1;
    end
    if length(attrlist)==1
        attrlistnew=[];
    end
end
if ~isempty(attrlistnew)
    attrlist=attrlistnew;
else
    attrlist=[];
end
D2{1}=D(D(:,splitattr)<=splitpoint,:);
D2{2}=D(D(:,splitattr)>splitpoint,:);
for j=1:2
    if isempty(D2{j})
        continue;
    else 
        CreateDecisionTree_Test(D2{j},attrlist,icount+1,iparenttag,j);
    end
end