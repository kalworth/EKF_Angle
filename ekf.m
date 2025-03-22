% 加载陀螺仪数据
% 顺序ax,ay,az,gx,gy,gz
% acc的单位是g
% gyro的单位是deg/s
data = load('D:/Serial Debug 2023-12-25 215322.csv');

% 获取数据个数 一组数据 ax ay az gx gy gz  
numRows = size(data, 1); 

% 迭代次数
epoch = 2;

% 欧拉角
global eul_deg;
eul_deg = zeros(numRows*epoch,3);

% 加速度计欧拉角
global eul_acc_deg;
global gz_yaw;
eul_acc_deg = zeros(numRows*epoch,3);
gz_yaw = 0;

% 更新周期
global period_T;
period_T = 0.001;

% 更新周期的一半
global Half_T;
Half_T = period_T / 2;

 % 机体四元数
global q; 
q = transpose([1 0 0 0]);

% 观测向量
global Z;

% 状态转移矩阵
global A;

% 观测矩阵
global H;

% 误差协方差矩阵
global P;
P_value = 100000;
P_diag = [P_value P_value P_value P_value];
P = diag(P_diag);

% 过程噪声协方差矩阵
global Q;
Q_value = 0.01;
Q_diag = [Q_value Q_value Q_value Q_value];
Q = diag(Q_diag);

% 测量噪声协方差矩阵
global R;
R_value = 10000000;
R_diag = [R_value R_value R_value];
R = diag(R_diag);

% 卡尔曼增益
global K;

global deg_num;
deg_num = 1;
for j = 1:epoch
    % 使用for循环遍历每一行  
    for i = 1:numRows  
        % 访问当前行的所有元素  
        row = data(i, :); % 使用冒号来选择所有列  

        ax = row(1) * 9.8;
        ay = row(2) * 9.8;
        az = row(3) * 9.8;
        gx = row(4) / 57.3;
        gy = row(5)/ 57.3;
        gz = row(6)/ 57.3;

        % 先验估计
        A =  update_A(gx,gy,gz);
        q = A * q; 
        q = q / norm(q);
        
        % 计算先验误差协方差
        AT = transpose(A);
        P = A*P*AT + Q;

        % 更新观测向量
        Z = transpose([ax,ay,az]);

        % 更新观测矩阵
        H = update_H(q);

        % 计算卡尔曼增益
        HT = transpose(H);
        K = (P*HT) / (H*P*HT+R);

        % 后验估计
        q = q + K * (Z - update_hx(q));
        q = q / norm(q);
        
        % 更新误差协方差
        P = (eye(4) - K * H) * P;

        % 欧拉角解算
        qt = transpose(q);
        eul = quat2eul(qt);  
        eul_deg(deg_num,:) = rad2deg(eul);

        gz_yaw = gz_yaw + gz*period_T;
        eul_acc_deg(deg_num,1) = rad2deg(atan2(ay,az));
        eul_acc_deg(deg_num,2) = -rad2deg(atan(ax/sqrt(ay*ay+az*az)));
        eul_acc_deg(deg_num,3) = gz_yaw;
        deg_num = deg_num + 1;
        
        norm(q);
    end
end

line_width = 2;
figure; % 激活第三个子图  
plot(1:numRows*epoch, eul_acc_deg(:,1), 'r-','LineWidth', line_width); % 绘制滚转角（roll）  
hold on; % 保持当前图形，以便添加更多曲线  
plot(1:numRows*epoch, eul_deg(:,3), 'b-','LineWidth', line_width); % 绘制滚转角（roll）  
legend('acc', 'ekf');  
title('Pitch Angle (degrees)');  
xlabel('Time');  
ylabel('Angle'); 

figure; % 激活第三个子图  
plot(1:numRows*epoch, eul_acc_deg(:,2), 'r-','LineWidth', line_width); % 绘制滚转角（roll）  
hold on; % 保持当前图形，以便添加更多曲线  
plot(1:numRows*epoch, eul_deg(:,2), 'b-','LineWidth', line_width); % 绘制滚转角（roll）  
legend('acc', 'ekf'); 
title('Roll Angle (degrees)');  
xlabel('Time');  
ylabel('Angle'); 


% 状态转移矩阵更新
function A_Update = update_A(gx,gy,gz)
    global Half_T; 
    gx_ = gx * Half_T;
    gy_ = gy * Half_T;
    gz_ = gz * Half_T;
    A_Update = [  1  -gx_ -gy_  -gz_;
                 gx_  1    gz_  -gy_;
                 gy_ -gz_   1    gx_; 
                 gz_  gy_ -gx_    1 ;];
end

% 测量矩阵更新
function H_Update = update_H(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    H_Update = [-2*q2  2*q3 -2*q0 2*q1;
                 2*q1  2*q0  2*q3 2*q2;
                 2*q0 -2*q1 -2*q2 2*q3];
end

% 观测加速度更新
function hx_Update = update_hx(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    hx_Update = [2*q1*q3-2*q0*q2;
                 2*q2*q3+2*q0*q1;
                 1-2*q1^2-2*q2^2];
end
