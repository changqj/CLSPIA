% 既插值又拟合  (迭代法2 我们的方法)
tic
clear
close all
%% 初始变量赋值
preserve_flag =0;
isPlot = 0;circle_r=0.22;%0.22;
kkk=0;
MaxStep = 100000; % 设置最大迭代步骤
mat_file = {'eg1','eg2','eg3',...
    'eg4','eg5','woody'};
interpolatedPoints = {[1:20:268],1:50:577,1:50:577,1:25:273,1:20:205,1:12:128,...
    1:20:205,1:10:97,1:20:180,[540,759,1191,1387,1494,1627,1967,2239,...
    2436,2512,2544,2552,2576,2859,3405,3501,3517,3805,4734,4985,5061,...
    5309,5473,5498,5679,5854,6162,6552,6992,7514,7539,7573,8144,8310,...
    8409,8910,9173,9294,9594,9599],[6 42 75 110 147 182],1:50:628,...
    [1 21 43 63 85 105 127 147 168 190],2*[10,20,36,58,69,79,93,110,...
    121,131,141,152,161,177,195,204,213,230,250],...
    [1 28 55 75 102 132 143 165 190 210 226 250 273 290 310 335 360]};% 每个例子对应的插值点列[20 28 244 260 275 280 291 326 340 564]

ctrl_point = [55,55,55,30,30,22,20,30,20,1000,37,200,22,50,60];
ii=3;% 例子 
%%%%%%%%%%% 插值点=控制点个数的例子 %%%%%%%%%%%%%%%%
te = [5 14 28 45 55 65 75 80 102 115 132 143 155 165 178 190 200 210 212 226 238 250 261 273 281 290 300 310 322 335 348 355 360 365 370];
interpolatedPoints{15} = te(unique([1:7:35 6 10 14 17 20 23 27 31 4 12 18 25 28 30 33]));
ctrl_point(15) = 35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chazhidian{ii} = 35:50:200;
load(mat_file{ii}) % 导入例子数据
dataPoints = dataPoints';
dataPoints1 = dataPoints;
R_index = interpolatedPoints{ii}; % 获得插值点列
Q_index = setdiff(1:size(dataPoints,2),R_index);  % 获得逼近点列
%%%%%%%%%%% 插值点=控制点个数的例子 %%%%%%%%%%%%%%%%
% dataPoints(:,Q_index) = dataPoints(:,Q_index) + 0.04*rands(1)*ones(size(dataPoints(:,Q_index)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = dataPoints(:,Q_index);  % 取出要拟合的数据点集
R = dataPoints(:,R_index);  % 取出要插值的数据点集
ParaOfDataPoint = Parameterize(dataPoints1',1);  % 将初始的几个数据点参数化到t
T = ParaOfDataPoint(Q_index);   % 得到拟合数据点对应的参数集
S = ParaOfDataPoint(R_index);   % 得到插值数据点对应的参数集

% 从数据点中选择控制顶点 LSPIA的方法
ctrl_point_num = ctrl_point(ii);  % 控制顶点个数% ctrl_point_num = floor(size(dataPoints,2)/10);
[P0,i_step,deviation,knot] = cqj_LSPIA(dataPoints1,ctrl_point_num,1000,1e-2,preserve_flag);  % 调用LSPIA方法计算初始控制顶点
% P0 = ones(size(P0));
lambda = zeros(size(R));      % 初始化插值点的lamda向量
%% 计算对称矩阵A^TA并得到miu
temp_A = eye(size(dataPoints,2),ctrl_point_num);
for i = 1:size(dataPoints,2)
    for j = 1:ctrl_point_num
        temp_A(i,j) = cqj_BsplineBasisFunction(j-1,3,ParaOfDataPoint(i),knot);
    end
end
A = temp_A(Q_index,:);
B = temp_A';
B = B(:,R_index);
ATA2 = 2*(A'*A);
ATQ2 = 2*A'*Q';
[~,eigenValue] = eig(B'/(ATA2)*B);
mu = 1/max(diag(eigenValue));
[~,eigenValueOfmiu1] = eig(ATA2);
mu1 = max(diag(eigenValueOfmiu1));


Error = [B'*P0'-R'];
E(1) = sum(Error(:).^2);
%% PLOT
if isPlot
    P=P0;
    figure('visible','off');
    box on;grid on;%title('既插值又拟合');
    axis tight;                     % 紧坐标轴
    axis equal;                     % 等比坐标轴
    axis([0 1 0 1]);
    hold on
    
    h1 = plot(dataPoints(1,R_index),dataPoints(2,R_index),'*b');    % 画插值点
    h2 = plot(dataPoints(1,Q_index),dataPoints(2,Q_index),'.b');    % 画逼近点
    
    r = [];
    for t = sort(unique([0:0.001:1,ParaOfDataPoint]))
        B1 = zeros(ctrl_point_num,1);
        for i = 0:ctrl_point_num-1
            B1(i+1) = cqj_BsplineBasisFunction(i,3,t,knot);
        end
        r =  [r P*B1];
    end
    %h3 = plot(P(1,:),P(2,:),'-k');                      % 画控制顶点
    h4 = plot(r(1,:),r(2,:),'-r');                      % 画3次B样条曲线
    % legend([h1 h2 h3 h4],'插值点','拟合点','控制多边形','3次B样条曲线');
    % legend([h1 h2 h4],'插值点','拟合点','3次B样条曲线');
    kkk=kkk+1;
    saveas(gcf,['res/' num2str(kkk) '_1.png']);
    hold off
end
%% 迭代
res = [];
itrNum = [];
for m = 1:MaxStep
    error0 = 0;
    for k = 0:10*MaxStep
        %     tempP= (temp1\(temp2-newB*lamda'))';
        errortemp = (1/mu1*(ATQ2-B*lambda'-ATA2*P0'))';
        P0 = P0+errortemp;
        tP0 = P0;
        %% 加保形约束
        P0 = preserve(P0,errortemp,preserve_flag);
        sum(sum(P0 ~= tP0));
        %%
        Error = [B'*P0'-R'];
        E(end+1) = sum(Error(:).^2);
        %% PLOT
        if isPlot
            P=P0;
            kkk=kkk+1;
            figure('visible','off');
            axis([0 1 0 1]);
            box on;grid on;%title('既插值又拟合');
            axis tight;                     % 紧坐标轴
            axis equal;                     % 等比坐标轴
            axis([0 1 0 1]);
            hold on
            
            h1 = plot(dataPoints(1,R_index),dataPoints(2,R_index),'*b');    % 画插值点
            h2 = plot(dataPoints(1,Q_index),dataPoints(2,Q_index),'.b');    % 画逼近点
            
            r = [];
            for t = sort(unique([0:0.001:1,ParaOfDataPoint]))
                B1 = zeros(ctrl_point_num,1);
                for i = 0:ctrl_point_num-1
                    B1(i+1) = cqj_BsplineBasisFunction(i,3,t,knot);
                end
                r =  [r P*B1];
            end
            
            %         h3 = plot(P(1,:),P(2,:),'-k');                      % 画控制顶点
            h4 = plot(r(1,:),r(2,:),'-r');                      % 画3次B样条曲线
            % legend([h1 h2 h3 h4],'插值点','拟合点','控制多边形','3次B样条曲线');
            % legend([h1 h2 h4],'插值点','拟合点','3次B样条曲线');
            saveas(gcf,['workPic/' num2str(kkk) '_1.png']);
            hold off
        end
        %%
        temp_1 = errortemp-error0;
        if max(diag(errortemp'*errortemp))<(1e-3)^2
            break;
        end
        error0 = errortemp;
    end
    itrNum(m) = k;
    %     if max(diag((P0-P)'*(P0-P)))<(1e-10)^2
    %         P = P0;
    %         break;
    %     end
    P = P0;
    res(:,:,end+1)=P;
    %     max(diag((B'*P'-R')*(B'*P'-R')'))
    if max(diag((B'*P'-R')*(B'*P'-R')'))<(1e-10)^2% && max(sum((A*P'-Q').^2,2))<(5e-2)^2
        break;
    end
    lambda = lambda+mu*(B'*P'-R')';
    
end

aError = [A*P0'-Q'];
        aE = sum(aError(:).^2);

TEMPA = [2*(A'*A) B;B' zeros(size(B,2))];
TEMPB = [2*A'*Q';R'];
TEMP = TEMPA\TEMPB;
P1 = TEMP(1:size(P,2),:);
P1 =P1';
% P0 = P1;
if ~isPlot
    figure
    %     magnify
    P=P0;
    box on;%grid on;%title('既插值又拟合');
    hold on
    axis tight;                     % 紧坐标轴
    axis equal;                     % 等比坐标轴
    %         axis([0 0.3 -0.1 0.1]);
    
    h1 = plot(dataPoints(1,R_index),dataPoints(2,R_index),'*r');    % 画插值点
    h2 = plot(dataPoints(1,Q_index),dataPoints(2,Q_index),'.b');    % 画逼近点
    
    r = [];
    for t = sort(unique([0:0.001:1,ParaOfDataPoint]))
        B1 = zeros(ctrl_point_num,1);
        for i = 0:ctrl_point_num-1
            B1(i+1) = cqj_BsplineBasisFunction(i,3,t,knot);
        end
        r =  [r P*B1];
    end
    
    %  h3 = plot(P(1,:),P(2,:),'-ok');                      % 画控制顶点
    h4 = plot(r(1,:),r(2,:),'-r');                      % 画3次B样条曲线
    % legend([h1 h2 h3 h4],'插值点','拟合点','控制多边形','3次B样条曲线');
    % legend([h1 h2 h4],'插值点','拟合点','3次B样条曲线');
    if 0
        O = dataPoints(:,R_index(9));
        %     theta = (0.25:0.5:2.25)*pi;
        theta = (0:0.01:2)*pi;
        rr = circle_r*ones(size(theta));
        X = rr.*cos(theta)+O(1);
        Y = rr.*sin(theta)+O(2);
        plot(X,Y,'-k')
        
        figure,hold on,axis tight;                     % 紧坐标轴
        axis equal;                     % 等比坐标轴
        axis off
        % part1
        part1 = [];
        for i = 1:size(R_index,2)
            if ~(norm(dataPoints(:,R_index(i))-O)>circle_r)
                part1(:,end+1) = dataPoints(:,R_index(i));
            end
        end
        % part2
        part2 = [];
        for i = 1:size(Q_index,2)
            if ~(norm(dataPoints(:,Q_index(i))-O)>circle_r)
                part2(:,end+1) = dataPoints(:,Q_index(i));
            end
        end
        % part4
        part4 = {};
        tempa = [];
        kkkk = 1;flag = 0;
        for i = 1:size(r,2)
            if ~(norm(r(:,i)-O)>circle_r)
                flag = 1;
                tempa(:,end+1) = r(:,i);
                if i==size(r,2)
                    part4{kkkk} = tempa;
                end
            elseif (norm(r(:,i)-O)>circle_r) && flag == 1
                part4{kkkk} = tempa;
                flag = 0;tempa = [];
                kkkk = kkkk + 1;
            end
        end
        part1 = part1;
        part2 = part2;
        
        plot(part1(1,:),part1(2,:),'*r')
        plot(part2(1,:),part2(2,:),'.b')
        for i = 1:size(part4,2)
            part4{i} = part4{i};
            plot(part4{i}(1,:),part4{i}(2,:),'-r');
        end
        plot(X,Y,'-k')
    end
end
toc