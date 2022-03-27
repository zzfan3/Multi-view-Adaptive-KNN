clear;clc;
  
%read different view data
[X,Y] = read_data_MSRC(1, 2);
class = 7; % in multikTreeclassify.m line 49 class also need to change
viewnum = 5;
[m1, n1] = size(X{1});
ind = crossvalind('Kfold', m1, 10);% ten-fold


% kTree code for "Efficient kNN Classification With Different Numbers of Nearest Neighbors"
for im = 1:10 
for i = 1:viewnum
test = (ind==im);    
train = ~test;
alltest{im} = test;

X_tr = X{i}(train,:);
X_te = X{i}(test,:);   
Y_tr = Y(train,1);
Y_te = Y(test,1);
trainnum = size(X_tr,1);  
testnum = size(X_te,1);  
[m, n] = size(X_tr);

tenum = randperm(size(X_te,1));
tennum1 = sort(tenum);
X_te = [tennum1', X_te];
%% prepare before decision tree
Xdt = [];Ydt = [];   
X_trconvert = X_tr';
for j = 1:trainnum
    Xdt{j} = X_tr';
    Ydt{j} = X_trconvert(:,j);
end
%% using training samples to reconstruct themselves
tic;
for ek = 1:3
    if ek==1
       rho1_ts=0;  rho4_ts=5*10^-1; rho3_ts=5*10^-5;
    elseif ek==2
        rho1_ts=0;  rho4_ts=5*10^-1; rho3_ts=5*10^-3;
    elseif ek==3
        rho1_ts=0;  rho4_ts=5*10^-3; rho3_ts=5*10^-1;
    end
p_ts = 1;
for j2_ts = 1:length(rho3_ts)
    for k2_ts = 1:length(rho4_ts)
        [W_ts, funcVal_ts] = L2LPP_L21L1_gpu(Xdt,Ydt,rho1_ts,rho3_ts(j2_ts),rho4_ts(k2_ts));
        W_ts(W_ts<0) = 0;    
        logicW_ts = W_ts&1;  
        colvalue_ts = sum(logicW_ts);
        if all(colvalue_ts)
            Wcollection_ts(p_ts) = {W_ts};
            indexCollection_ts(p_ts) = {[j2_ts k2_ts]};
            p_ts = p_ts+1;
        end
    end
end
%% build kTree
global icount;    
global attrNodeDC;   
attrNodeDC = struct([]);
icount = 0;
attrlist = 1:1:n;
pn_ts = 1;
w_ts = Wcollection_ts{pn_ts};         %w:n*c
logicw_ts = w_ts&1; 
colvalue_ts = sum(logicw_ts);
YDs = colvalue_ts';
D = [X_tr, YDs];
AttrNodeDC0 = CreateDecisionTree_Test(D,attrlist,1,1,0);        % build decision tree  
AttrNodeDC = FillMissingChildNode_DC20151021(1,AttrNodeDC0);      % Supplement the missing child node
%% use kTree to classify
% tic;
Xtrain = [X_tr, Y_tr];
Xtest = [X_te, Y_te];
[row, col] = size(Xtrain);
classnum = class;
predict = multikTreeClassify(AttrNodeDC,1,Xtrain,row,col,Xtest,classnum); % use kTree to classify test samples
predict = sortrows(predict,1);

predictek{ek} = predict;
for jt = 1:testnum
    Eek{ek}(jt,:) = predict(jt,2:1+class)/sum(predict(jt,2:1+class));
    [~, v_ek{ek}(jt)] = max(Eek{ek}(jt,:));
    vek_predict{i,im}{ek}(jt) = v_ek{ek}(jt);
end
vek_INDF = find((vek_predict{i,im}{ek}'-Y_te)==0);
vek_correctnum(ek) = length(vek_INDF)/testnum*100
end

% multi_view evidence 2022.3
[~, mek] = max(vek_correctnum);
predict = predictek{mek};
% calculate evidence parameters
for jt = 1:testnum
    % D-S combine
    E{i,im}(jt,:) = predict(jt,2:1+class)/sum(predict(jt,2:1+class));
    alpha{i,im}(jt,:) = E{i,im}(jt,:)+1;
    S{i,im}(jt) = sum(alpha{i,im}(jt,:));
    b{i,im}(jt,:) = E{i,im}(jt,:)/S{i,im}(jt);
    u{i,im}(jt) = class/S{i,im}(jt);
    
    [~, v_e{i,im}(jt)] = max(E{i,im}(jt,:));
    v_predict{i,im}(jt) = v_e{i,im}(jt);
end

% calculate correctnum
v_INDF=find((v_predict{i,im}'-Y_te)==0);
v_correctnum=length(v_INDF);

% im_th interation, i_th view
correct_vkTree(i,im)=v_correctnum/testnum*100
rmse_kTree(i,im)=rmse(Y_te,v_predict{i,im}')
end

% combine evidence
[multi_predict]=CombineEvidence(b,u,testnum,viewnum,im,class);

% calculate the correct number
INDF=find((multi_predict{im}'-Y_te)==0);
multi_correctnum=length(INDF);
correct_multikTree(im)=multi_correctnum/testnum*100  

end

mcorrect_vkTree=mean(correct_vkTree,2)'% each view's mean correctnum

% all view's max correctnum
all_mcorrect_vkTree=max(mcorrect_vkTree)
mcorrect_multikTree = mean(correct_multikTree)











