clear;clc;

A=load('German.txt');  %    download dataset
A=A(1:250,:);
XdataNor=A(:,1:24);    %    training samples without labels
XdataNor= NormalizeFea(XdataNor,0);
[m n]=size(XdataNor);
cmaker=A(:,25);        %    labels

ind=crossvalind('Kfold',m,10);

for i=1:10
test=(ind==i);    
train=~test;

alltest{i}=test;

X_tr=XdataNor(train,:);   
X_te=XdataNor(test,:);   

Y_tr=cmaker(train,1);
Y_te=cmaker(test,1);
trainnum=size(X_tr,1);  
testnum=size(X_te,1);  
%% prepare before decision tree
trainnum=size(X_tr,1);
Xdt=[];Ydt=[];   
X_trconvert=X_tr';
for j=1:trainnum
    Xdt{j}=X_tr';
    Ydt{j}=X_trconvert(:,j);
end

%% using training samples to reconstruct themselves
% German dataset parameter settings
tic;
rho1_ts=0;  rho4_ts=5*10^-3; rho3_ts=7.5*10^-3;    % L21         2015.10.26
p_ts=1;
for j2_ts=1:length(rho3_ts)
    for k2_ts=1:length(rho4_ts)
        [W_ts, funcVal_ts] = L2LPP_L21L1(Xdt,Ydt,rho1_ts,rho3_ts(j2_ts),rho4_ts(k2_ts));
        W_ts(W_ts<0)=0;    
        logicW_ts=W_ts&1;  
        colvalue_ts=sum(logicW_ts);
        if all(colvalue_ts)
            Wcollection_ts(p_ts)={W_ts};
            indexCollection_ts(p_ts)={[j2_ts k2_ts]};
            p_ts=p_ts+1;
        end
    end
end
%% build kTree
global icount;    
global attrNodeDC;   
attrNodeDC=struct([]);
icount=0;
attrlist=1:1:n;
pn_ts=1;
w_ts=Wcollection_ts{pn_ts};         %w:n*c
logicw_ts=w_ts&1; 
colvalue_ts=sum(logicw_ts);
YDs=colvalue_ts';
D=[X_tr YDs];
AttrNodeDC0=CreateDecisionTree_Test(D,attrlist,1,1,0);        % build decision tree  
AttrNodeDC=FillMissingChildNode_DC20151021(1,AttrNodeDC0);      % Supplement the missing child node
%% use kTree to classify
tic;
Xtrain=[X_tr Y_tr];
Xtest=[X_te Y_te];
[row col]=size(Xtrain);
ClassfiyResultkTree=kTreeClassify20151021(AttrNodeDC,1,Xtrain,row,col,Xtest); % use kTree to classify test samples
correctnumkTree=ClassfiyResultkTree;
correct_kTree(i)=correctnumkTree/testnum;                 % classification accuracy
time_kTree(i)=toc;                                        % classification time
end
mcorrect_kTree=mean(correct_kTree);

