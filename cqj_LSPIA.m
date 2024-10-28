function [P,i_step,meanDeviation,knot] = cqj_LSPIA(init_data,...
    ctrl_point_num,max_step,accuracy,preserve_flag)
%{
% 最小二乘迭代法（LSPIA）
% init_data : 初始控制点，如果是三维点，size(data) = [3 n]
% ctrl_point_num：控制点的个数
% max_step : 最大迭代次数
% accuracy : 精度
% P : 最后一次的控制点
% i_step : 迭代次数
% deviation : 结果曲线上的对应点与初始控制点的偏差
%}

if exist('preserve_flag', 'var') ~= 1
    preserve_flag = 0;
end
if ctrl_point_num > size(init_data,2)
    error('控制顶点数量过大！\n');
end
meanDeviation = 0;
Q = init_data;

%% 从数据点中选择控制顶点
P(:,1) = Q(:,1);
P(:,ctrl_point_num) = Q(:,end);
m = size(init_data,2) - 1;
n = ctrl_point_num - 1;
for i = 2:n
    P(:,i) = Q(:,floor((m+1)*(i-1)/n)+1);
end
%% 参数化并得到节点矢量
ParaOfDataPoint = Parameterize(Q',1);   % 将初始的几个数据点参数化到t
d = (m+1)/(n-2);
knot = zeros(1,n+1);
for j = 1:n-3
    i = floor(j*d);
    alpha = j*d-i;
    knot(j+4) = (1-alpha)*ParaOfDataPoint(i)+alpha*ParaOfDataPoint(i+1);
end
knot = [knot ones(1,4)];

%% 计算对称矩阵A^TA并得到miu
A = zeros(m+1,n+1);
for i = 1:m+1
    for j = 1:n+1
        A(i,j) = cqj_BsplineBasisFunction(j-1,3,ParaOfDataPoint(i),knot);
    end
end
miu = 2/max(sum(A'*A,2));

%% 循环迭代
for i_step = 1:max_step
    % 计算小derta
    P_ = P*A';
    derta = Q - P_;
    % 用数据点与曲线上点的距离判断是否到达精度
    deviation = vecnorm(derta);
    meanDeviation = (mean(deviation))^0.5;
    Derta = miu*(derta*A);   % 计算大derta
%     hold off
%     plot(P(1,:),P(2,:),'-or');
%     hold on
    P = P + Derta;  % 更新控制顶点
%     plot(P(1,:),P(2,:),'-xg');
    %% 加入保形约束
    P = preserve(P,Derta,preserve_flag);
    %%
%     plot(P(1,:),P(2,:),'-vb');hold on
%     plot(P(1,[8 10]),P(2,[8 10]),'--');
    if mean(diag(Derta'*Derta))^0.5 < accuracy     % 如果精度满足则跳出循环
        break;
    end
end
end