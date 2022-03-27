function [X,Y]=read_data_MSRC(epoch, t)

data = load('./data/MSRC-v1.mat');

X1 = data.X{1,1};
X2 = data.X{1,2};
X3 = data.X{1,3};
X4 = data.X{1,4};
X5 = data.X{1,5};
Y = data.Y;

X{1} = X1(epoch:t:end,:);
X{2} = X2(epoch:t:end,:);
X{3} = X3(epoch:t:end,:);
X{4} = X4(epoch:t:end,:);
X{5} = X5(epoch:t:end,:);
Y = Y(epoch:t:end,:);

