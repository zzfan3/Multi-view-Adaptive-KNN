%% Classify multiple samples 
%% input: Xtrain denotes X_tr+Y_tr; Xtest denotes X_te+Y_te
%% output:predictions
function [mpredictY]=multikTreeClassify(attrNodekTree,i,Xtrain,row,col,Xtest,Xtestrow,class)

if ~isempty(attrNodekTree(1,i).Lleaflabel)
    kval=attrNodekTree(1,i).Lleaflabel;
    [rr, rc]=size(Xtest); % Xtest第一列是序号
    for it=1:rr
        if Xtest(it,rc)~=0
            ind_kNN=knnsearch(Xtrain(:,1:col-1),Xtest(it,2:rc-1),'k',kval);%第it个测试样本的k近邻
            sample=Xtrain(ind_kNN,:);
            samplelabel=sample(:,col);
            % calculate e
            c(it,:)=zeros(1,class);
            for i=1:class
                c(it,i)=length(find(samplelabel==i));
            end
            predictY(it,:)=[Xtest(it,1),c(it,:)]; % smaple number and Knn labels
        end
    end
    mpredictY=predictY;    
    return;
end
if  ~isempty(attrNodekTree(1,i).Rleaflabel)
    kval=attrNodekTree(1,i).Rleaflabel;
    [rr, rc]=size(Xtest);
    for it=1:rr
        if Xtest(it,rc)~=0
            ind_kNN=knnsearch(Xtrain(:,1:col-1),Xtest(it,2:rc-1),'k',kval);%第it个测试样本的k近邻
            sample=Xtrain(ind_kNN,:);
            samplelabel=sample(:,col);
            % calculate e
            c(it,:)=zeros(1,class);
            for i=1:class
                c(it,i)=length(find(samplelabel==i));
            end
            predictY(it,:)=[Xtest(it,1),c(it,:)]; % smaple number and Knn labels
        end
    end
    mpredictY=predictY; 
    return;
end

D{1}=Xtest(Xtest(:,1+attrNodekTree(1,i).splitattr)<=attrNodekTree(1,i).splitpoint,:);%Xtest多了一列
D{2}=Xtest(Xtest(:,1+attrNodekTree(1,i).splitattr)>attrNodekTree(1,i).splitpoint,:);

% need to change
class=7;
%-------------

for j=1:2
    if isempty(D{j}) %all(D{j})
        if j==1
            LLpredictlabel=[];
        end
        if j==2
            RRpredictlabel=[];
        end
        continue;
    end
    if j==1
        LLpredictlabel=multikTreeClassify(attrNodekTree,attrNodekTree(1,i).leftchildNode,Xtrain,row,col,D{1},Xtestrow,class);
        LLpredictlabel=double(LLpredictlabel);
    end
    if j==2
        RRpredictlabel=multikTreeClassify(attrNodekTree,attrNodekTree(1,i).rightchildNode,Xtrain,row,col,D{2},Xtestrow,class);
        RRpredictlabel=double(RRpredictlabel);
    end
end
mpredictY=vertcat(LLpredictlabel,RRpredictlabel);
