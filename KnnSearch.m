function [index]=KnnSearch(Xtr,xte,k)
[m n]=size(Xtr);
aa=sum(Xtr.*Xtr,2);
bb=sum(xte.*xte,2);
ab=Xtr*xte';
D2=aa+bb-2.*ab;
D=sqrt(D2);
SD=sort(D);
for i=1:k
    indtemp=find(SD(i)==D)
    index(i)=indtemp(1);
end