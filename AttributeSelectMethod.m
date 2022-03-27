%% function: Attribute selection method (information gain): for continuous value attributes
%% input: D denotes training samples + labels, attrlist is used for record weather the atrribute is used or not, 
%%        initial values are [1 2 3 4 ....], denotes [attribute1, attribute2, attrbute3,...]
%% output:splitattr denotes the selected split atrribute, and splitpoint denotes split value (optimal splitpoint)
function [splitattr, splitpoint]=AttributeSelectMethod(D,attrlist)
%
[r c]=size(D);
X=D(:,1:c-1);
Y=D(:,c);
cateY=unique(Y);
freqY=histc(Y,cateY);
rownum=size(X,1);
for i=1:length(attrlist)
    attr=attrlist(i);
    [sortattr sortind]=sort(X(:,attr));
    for j=1:rownum-1
        maysplitpoint(j)=(sortattr(j)+sortattr(j+1))/2;
        
        leftY=Y(sortind(1:j,:));
        catelY=unique(leftY);
        freqlY=histc(leftY,catelY);
        pl=freqlY./j;
        suml=0;
        for kl=1:length(catelY)
            suml=-pl(kl)*log2(pl(kl))+suml;
        end
        rightY=Y(sortind(j+1:rownum,:));
        caterY=unique(rightY);
        freqrY=histc(rightY,caterY);
        pr=freqrY./(rownum-j);
        sumr=0;
        for kr=1:length(caterY)
            sumr=-pr(kr)*log2(pr(kr))+sumr;
        end
        maysplitpointInfor(j)=j/rownum*suml+(rownum-j)/rownum*sumr;
    end
    [minInforval(i) minInforind(i)]=min(maysplitpointInfor);
    eachsplitpoint(i)=maysplitpoint(minInforind(i));
end
[minInfor minind]=min(minInforval);
splitattr=attrlist(minind);
splitpoint=eachsplitpoint(minind);