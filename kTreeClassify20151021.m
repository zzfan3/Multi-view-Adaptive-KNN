%% Classify multiple samples 
%% input: Xtrain denotes X_tr+Y_tr; Xtest denotes X_te+Y_te
%% output:The correct number of predictions
function [correctnum]=kTreeClassify20151021(attrNodekTree,i,Xtrain,row,col,Xtest)

if ~isempty(attrNodekTree(1,i).Lleaflabel)
    kval=attrNodekTree(1,i).Lleaflabel;
    [rr rc]=size(Xtest);
    for it=1:rr
        ind_kNN=knnsearch(Xtrain(:,1:col-1),Xtest(it,1:rc-1),'k',kval);%第it个测试样本的k近邻
        sample=Xtrain(ind_kNN,:);
        samplelabel=sample(:,col);
        ul=unique(samplelabel);
        fy=histc(samplelabel,ul);
        [mv mi]=max(fy);
        predictlabel(it)=ul(mi);
    end
    INDF=find((predictlabel'-Xtest(:,rc))==0);
    correctnum=length(INDF); 
    return;
end
if  ~isempty(attrNodekTree(1,i).Rleaflabel)
    kval=attrNodekTree(1,i).Rleaflabel;
    [rr rc]=size(Xtest);
    for it=1:rr
        ind_kNN=knnsearch(Xtrain(:,1:col-1),Xtest(it,1:rc-1),'k',kval);%第it个测试样本的k近邻
        sample=Xtrain(ind_kNN,:);
        samplelabel=sample(:,col);
        ul=unique(samplelabel);
        fy=histc(samplelabel,ul);
        [mv mi]=max(fy);
        predictlabel(it)=ul(mi);
    end
    INDF=find((predictlabel'-Xtest(:,rc))==0);
    correctnum=length(INDF); 
    return;
end
D{1}=Xtest(Xtest(:,attrNodekTree(1,i).splitattr)<=attrNodekTree(1,i).splitpoint,:);
D{2}=Xtest(Xtest(:,attrNodekTree(1,i).splitattr)>attrNodekTree(1,i).splitpoint,:);
for j=1:2
    if isempty(D{j})
        if j==1
            Lcorrectnum=0;
        end
        if j==2
            Rcorrectnum=0;
        end
        continue;
    end
    if j==1
        Lcorrectnum=kTreeClassify20151021(attrNodekTree,attrNodekTree(1,i).leftchildNode,Xtrain,row,col,D{j});
    end
    if j==2
        Rcorrectnum=kTreeClassify20151021(attrNodekTree,attrNodekTree(1,i).rightchildNode,Xtrain,row,col,D{j});
    end
end
correctnum=Lcorrectnum+Rcorrectnum;