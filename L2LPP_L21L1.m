%0.5 * norm (Y{1,1} - X{ 1,1}* W(:,1))^2
% 0.5*rho3 * W(:,1)'*L{1,1}*W(:,1)
function [W, funcVal] = L2LPP_L21L1(X,Y,rho1,rho3,rho4,C,opts)

%% OBJECTIVE
% argmin_W {sum_i^t (0.5 * norm (Y{i} - X{i} * W(:, i))^2)
%            + opts.rho_L2 * \|W\|_2^2 + rho1 * \|W\|_{2,1} }
%            + rho3* tr(W'*L*W) + rho4*\|W\|_1    % 修改为 tr(W'X'LXW)
%% INPUT
% X: {d * n} * t - input matrix
% C: C{i}=X{i}'*L*X{i}, Laplacian matrix
% rho1: L2,1-norm group Lasso parameter.
% rho4: L1-norm group Lasso parameter.
% optional:
%   opts.rho_L2: L2-norm parameter (default = 0).
%   opts.rho_L3: manifold regularization parameter.
%
%% OUTPUT
% W: model: n * t or n * n
% funcVal: function value vector.

%% ------------------------ NOTE THAT ------------------------
%   Last modified on Nov. 11, 2013. 
%% ------------------------------------------------------------------------

%% example for reconstruction itself
%clear;clc;XX = rand(252,700);for i=1:700,X{i} = XX;end;for i=1:700,Y{i} = XX(:,i);end;
%[W, funcVal] = L2LPP_L21L1(X, Y,1,1,2);
%% example for reconstruction others, such as reconstructing test data Y with training data X
%clear;clc;XX = rand(252,700);for i=1:16,X{i} = XX;end;for i=1:16,Y{i} = rand(252,1);end;
%[W, funcVal] = L2LPP_L21L1(X, Y,1,1,2);

%genpath('C:\Users\xiaofeng\nips2012release\MALSAR1.1\MALSAR\')
%% Code starts here
if (~exist('C','var'))
    options = [];
    options.Metric = 'Euclidean';
    options.NeighborMode = 'KNN';
    options.k = 3;
    options.WeightMode = 'HeatKernel';
    options.t = 1;
    for ii=1:length(X)
        L{ii} = constructW(X{ii},options); % 修改过，原来代码： L{ii} = constructW(X{ii}',options);
%         C{ii} = L{ii};
        C{ii}=X{ii}'* L{ii}*X{ii}; % 修改过，原来没有这一行代码 
    end
end
% for ii=1:length(X)
%     C{ii} = L{ii};
% %     C{ii}=X{ii}'* L{ii}*X{ii}; % 修改过，原来没有这一行代码 
% end
if (~exist('opts','var'))
    opts.init = 0;       % guess start point from data.
    opts.tFlag = 1;      % terminate after relative objective value does not changes much.
    opts.tol = 10^-5;    % tolerance.
    opts.maxIter = 1000; % maximum iteration number of optimization.
end

if nargin <3
    error('\n Inputs: X, Y, rho1, should be specified!\n');
end
%X = multi_transpose(X);

if nargin <4
    opts = [];
end

% initialize options.
opts=init_opts(opts);

if isfield(opts, 'rho_L2')
    rho_L2 = opts.rho_L2;
else
    rho_L2 = 0;
end

% parameter for the graph regularizer
if (~exist('rho3','var'))
    rho3 = 1;
end

if (~exist('rho4','var'))
    rho4 = 2;
end


task_num  = length (X);
dimension = size(X{1}, 1);
funcVal = [];

XY = cell(task_num, 1);
W0_prep = [];
for t_idx = 1: task_num
    XY{t_idx} = X{t_idx}'*Y{t_idx};
    W0_prep = cat(2, W0_prep, XY{t_idx});
end

% initialize a starting point
if opts.init==2
    W0 = zeros(dimension, task_num);
elseif opts.init == 0
    W0 = W0_prep;
else
    if isfield(opts,'W0')
        W0=opts.W0;
        if (nnz(size(W0)-[dimension, task_num]))
            error('\n Check the input .W0');
        end
    else
        W0=W0_prep;
    end
end

bFlag=0; % this flag tests whether the gradient step only changes a little


Wz= W0;
Wz_old = W0;

t = 1;
t_old = 0;

iter = 0;
gamma = 1;
gamma_inc = 2;

while iter < opts.maxIter
    alpha = (t_old - 1) /t;
    
    Ws = (1 + alpha) * Wz - alpha * Wz_old;
    
    % compute function value and gradients of the search point
    gWs  = gradVal_eval(Ws);
    Fs   = funVal_eval (Ws);
    
    while true
        
        
        if rho4 == 0
            Wzp = FGLasso_projection(Ws - gWs/gamma, rho1 / gamma);
        else
            Us= l1_projection(Ws - gWs/gamma, 2 *  rho4 / gamma);
            Wzp = FGLasso_projection(Us, rho1 / gamma);
        end
       
        
        Fzp = funVal_eval  (Wzp);
        
        delta_Wzp = Wzp - Ws;
        r_sum = norm(delta_Wzp, 'fro')^2;
        Fzp_gamma = Fs + sum(sum(delta_Wzp.* gWs))...
            + gamma/2 * norm(delta_Wzp, 'fro')^2;
        
        if (r_sum <=1e-20)
            bFlag=1; % this shows that, the gradient step makes little improvement
            break;
        end
        
        if (Fzp <= Fzp_gamma)
            break;
        else
            gamma = gamma * gamma_inc;
        end
    end
    
    Wz_old = Wz;
    Wz = Wzp;
    
    funcVal = cat(1, funcVal, Fzp + nonsmooth_eval(Wz, rho1));
    
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end
    
    % test stop condition.
    switch(opts.tFlag)
        case 0
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iter>=2
                if (abs( funcVal(end) - funcVal(end-1) ) <=...
                        opts.tol* funcVal(end-1))
                    break;
                end
            end
        case 2
            if ( funcVal(end)<= opts.tol)
                break;
            end
        case 3
            if iter>=opts.maxIter
                break;
            end
    end
    
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
    
end

W = Wzp;

% private functions

    function [Wp] = FGLasso_projection (W, lambda )
        
        Wp = zeros(size(W));
        
        if opts.pFlag
            parfor i = 1 : size(W, 1)
                v = W(i, :);
                nm = norm(v, 2);
                if nm == 0
                    w = zeros(size(v));
                else
                    w = max(nm - lambda, 0)/nm * v;
                end
                Wp(i, :) = w';
            end
        else
            for i = 1 : size(W, 1)
                v = W(i, :);
                nm = norm(v, 2);
                if nm == 0
                    w = zeros(size(v));
                else
                    w = max(nm - lambda, 0)/nm * v;
                end
                Wp(i, :) = w';
            end
        end
    end


    function [z] = l1_projection (v, beta)
        
        z = zeros(size(v));
        vp = v - beta/2;
        z (v> beta/2)  = vp(v> beta/2);
        vn = v + beta/2;
        z (v< -beta/2) = vn(v< -beta/2);
        
        %l1_comp_val = sum(sum(abs(z)));
    end

    function [grad_W] = gradVal_eval(W)
        if opts.pFlag
            grad_W = zeros(zeros(W));
            parfor i = 1:task_num
                grad_W (i, :) = X{i}*(X{i}' * W(:,i)-Y{i});
            end
        else
            grad_W = [];
            for i = 1:task_num                
                grad_W = cat(2, grad_W, X{i}'*(X{i} * W(:,i)-Y{i})+ rho3*C{i}*W(:,i) );                
            end
        end
        grad_W = grad_W+ rho_L2 * 2 * W;
    end

% smooth part function value.
    function [funcVal] = funVal_eval (W)
        funcVal = 0;
        if opts.pFlag
            parfor i = 1: task_num                
                funcVal = funcVal + 0.5 * norm (Y{i} - X{i} * W(:, i))^2+ 0.5*rho3 * W(:,i)'*C{i}*W(:,i);                
            end
        else
            for i = 1: task_num                
                funcVal = funcVal + 0.5 * norm (Y{i} - X{i} * W(:, i))^2 + 0.5*rho3 * W(:,i)'*C{i}*W(:,i);               
            end
        end
        funcVal = funcVal + rho_L2 * norm(W,'fro')^2;
    end

    function [non_smooth_value] = nonsmooth_eval(W, rho_1)
        non_smooth_value = 0;
        if opts.pFlag
            parfor i = 1 : size(W, 1)
                w = W(i, :);
                non_smooth_value = non_smooth_value ...
                    + rho_1 * norm(w, 2);
            end
        else
            for i = 1 : size(W, 1)
                w = W(i, :);
                non_smooth_value = non_smooth_value ...
                    + rho_1 * norm(w, 2);
            end
        end        
        non_smooth_value=non_smooth_value+rho4*sum(sum(abs(W)));
    end

    function X = multi_transpose (X)
        for i = 1:length(X)
            X{i} = X{i}';
        end
    end

function opts = init_opts (opts)

%% Default values
DEFAULT_MAX_ITERATION = 1000;
DEFAULT_TOLERANCE     = 1e-4;
MINIMUM_TOLERANCE     = eps * 100;
DEFAULT_TERMINATION_COND = 1;
DEFAULT_INIT = 2;
DEFAULT_PARALLEL_SWITCH = false;

%% Starting Point
if isfield(opts,'init')
    if (opts.init~=0) && (opts.init~=1) && (opts.init~=2)
        opts.init=DEFAULT_INIT; % if .init is not 0, 1, or 2, then use the default 0
    end
    
    if (~isfield(opts,'W0')) && (opts.init==1)
        opts.init=DEFAULT_INIT; % if .W0 is not defined and .init=1, set .init=0
    end
else
    opts.init = DEFAULT_INIT; % if .init is not specified, use "0"
end

%% Tolerance
if isfield(opts, 'tol')
    % detect if the tolerance is smaller than minimum
    % tolerance allowed.
    if (opts.tol <MINIMUM_TOLERANCE)
        opts.tol = MINIMUM_TOLERANCE;
    end
else
    opts.tol = DEFAULT_TOLERANCE;
end

%% Maximum iteration steps
if isfield(opts, 'maxIter')
    if (opts.maxIter<1)
        opts.maxIter = DEFAULT_MAX_ITERATION;
    end
else
    opts.maxIter = DEFAULT_MAX_ITERATION;
end

%% Termination condition
if isfield(opts,'tFlag')
    if opts.tFlag<0
        opts.tFlag=0;
    elseif opts.tFlag>3
        opts.tFlag=3;
    else
        opts.tFlag=floor(opts.tFlag);
    end
else
    opts.tFlag = DEFAULT_TERMINATION_COND;
end

%% Parallel Options
if isfield(opts, 'pFlag')
    if opts.pFlag == true && ~exist('matlabpool', 'file')
        opts.pFlag = false;
        warning('MALSAR:PARALLEL','Parallel Toolbox is not detected, MALSAR is forced to turn off pFlag.');
    elseif opts.pFlag ~= true && opts.pFlag ~= false
        % validate the pFlag.
        opts.pFlag = DEFAULT_PARALLEL_SWITCH;
    end
else
    % if not set the pFlag to default.
    opts.pFlag = DEFAULT_PARALLEL_SWITCH;
end

if opts.pFlag
    % if the pFlag is checked,
    % check segmentation number.
    if isfield(opts, 'pSeg_num')
        if opts.pSeg_num < 0
            opts.pSeg_num = matlabpool('size');
        else
            opts.pSeg_num = ceil(opts.pSeg_num);
        end
    else
        opts.pSeg_num = matlabpool('size');
    end
end
end

end