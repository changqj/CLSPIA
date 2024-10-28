function [P,i_step,meanDeviation,knot] = cqj_LSPIA(init_data,...
    ctrl_point_num,max_step,accuracy,preserve_flag)
%{
% ��С���˵�������LSPIA��
% init_data : ��ʼ���Ƶ㣬�������ά�㣬size(data) = [3 n]
% ctrl_point_num�����Ƶ�ĸ���
% max_step : ����������
% accuracy : ����
% P : ���һ�εĿ��Ƶ�
% i_step : ��������
% deviation : ��������ϵĶ�Ӧ�����ʼ���Ƶ��ƫ��
%}

if exist('preserve_flag', 'var') ~= 1
    preserve_flag = 0;
end
if ctrl_point_num > size(init_data,2)
    error('���ƶ�����������\n');
end
meanDeviation = 0;
Q = init_data;

%% �����ݵ���ѡ����ƶ���
P(:,1) = Q(:,1);
P(:,ctrl_point_num) = Q(:,end);
m = size(init_data,2) - 1;
n = ctrl_point_num - 1;
for i = 2:n
    P(:,i) = Q(:,floor((m+1)*(i-1)/n)+1);
end
%% ���������õ��ڵ�ʸ��
ParaOfDataPoint = Parameterize(Q',1);   % ����ʼ�ļ������ݵ��������t
d = (m+1)/(n-2);
knot = zeros(1,n+1);
for j = 1:n-3
    i = floor(j*d);
    alpha = j*d-i;
    knot(j+4) = (1-alpha)*ParaOfDataPoint(i)+alpha*ParaOfDataPoint(i+1);
end
knot = [knot ones(1,4)];

%% ����Գƾ���A^TA���õ�miu
A = zeros(m+1,n+1);
for i = 1:m+1
    for j = 1:n+1
        A(i,j) = cqj_BsplineBasisFunction(j-1,3,ParaOfDataPoint(i),knot);
    end
end
miu = 2/max(sum(A'*A,2));

%% ѭ������
for i_step = 1:max_step
    % ����Сderta
    P_ = P*A';
    derta = Q - P_;
    % �����ݵ��������ϵ�ľ����ж��Ƿ񵽴ﾫ��
    deviation = vecnorm(derta);
    meanDeviation = (mean(deviation))^0.5;
    Derta = miu*(derta*A);   % �����derta
%     hold off
%     plot(P(1,:),P(2,:),'-or');
%     hold on
    P = P + Derta;  % ���¿��ƶ���
%     plot(P(1,:),P(2,:),'-xg');
    %% ���뱣��Լ��
    P = preserve(P,Derta,preserve_flag);
    %%
%     plot(P(1,:),P(2,:),'-vb');hold on
%     plot(P(1,[8 10]),P(2,[8 10]),'--');
    if mean(diag(Derta'*Derta))^0.5 < accuracy     % �����������������ѭ��
        break;
    end
end
end