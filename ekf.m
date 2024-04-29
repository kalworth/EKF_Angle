% ��������������
% ˳��ax,ay,az,gx,gy,gz
% acc�ĵ�λ��g
% gyro�ĵ�λ��deg/s
data = load('D:/Serial Debug 2023-12-25 215322.csv');

% ��ȡ���ݸ��� һ������ ax ay az gx gy gz  
numRows = size(data, 1); 

% ��������
epoch = 2;

% ŷ����
global eul_deg;
eul_deg = zeros(numRows*epoch,3);

% ���ٶȼ�ŷ����
global eul_acc_deg;
eul_acc_deg = zeros(numRows*epoch,2);

% ��������
global period_T;
period_T = 0.001;

% �������ڵ�һ��
global Half_T;
Half_T = period_T / 2;

 % ������Ԫ��
global q; 
q = transpose([1 0 0 0]);

% �۲�����
global Z;

% ״̬ת�ƾ���
global A;

% �۲����
global H;

% ���Э�������
global P;
P_value = 100000;
P_diag = [P_value P_value P_value P_value];
P = diag(P_diag);

% ��������Э�������
global Q;
Q_value = 0.01;
Q_diag = [Q_value Q_value Q_value Q_value];
Q = diag(Q_diag);

% ��������Э�������
global R;
R_value = 1000000;
R_diag = [R_value R_value R_value];
R = diag(R_diag);

% ����������
global K;

global deg_num;
deg_num = 1;
for j = 1:epoch
    % ʹ��forѭ������ÿһ��  
    for i = 1:numRows  
        % ���ʵ�ǰ�е�����Ԫ��  
        row = data(i, :); % ʹ��ð����ѡ��������  

        ax = row(1) * 9.8;
        ay = row(2) * 9.8;
        az = row(3) * 9.8;
        gx = row(4) / 57.3;
        gy = row(5)/ 57.3;
        gz = row(6)/ 57.3;

        % �������
        A =  update_A(gx,gy,gz);
        q = A * q; 

        % �����������Э����
        AT = transpose(A);
        P = A*P*AT + Q;

        % ���¹۲�����
        Z = transpose([ax,ay,az]);

        % ���¹۲����
        H = update_H(q);

        % ���㿨��������
        HT = transpose(H);
        K = (P*HT) / (H*P*HT+R);

        % �������
        q = q + K * (Z - H*q);

        % �������Э����
        P = (eye(4) - K * H) * P;

        % ŷ���ǽ���
        qt = transpose(q);
        eul = quat2eul(qt);  
        eul_deg(deg_num,:) = rad2deg(eul);

        eul_acc_deg(deg_num,1) = rad2deg(atan2(ay,az));
        eul_acc_deg(deg_num,2) = -rad2deg(atan(ax/sqrt(ay*ay+az*az)));
        deg_num = deg_num + 1;
    end
end

line_width = 2;
figure; % �����������ͼ  
plot(1:numRows*epoch, eul_deg(:,3), 'b-','LineWidth', line_width); % ���ƹ�ת�ǣ�roll��  
hold on; % ���ֵ�ǰͼ�Σ��Ա���Ӹ�������  
plot(1:numRows*epoch, eul_acc_deg(:,1), 'r-','LineWidth', line_width); % ���ƹ�ת�ǣ�roll��  
legend('ekf', 'acc');  
title('Pitch Angle (degrees)');  
xlabel('Time');  
ylabel('Angle'); 

figure; % �����������ͼ  
plot(1:numRows*epoch, eul_deg(:,2), 'b-','LineWidth', line_width); % ���ƹ�ת�ǣ�roll��  
hold on; % ���ֵ�ǰͼ�Σ��Ա���Ӹ�������  
plot(1:numRows*epoch, eul_acc_deg(:,2), 'r-','LineWidth', line_width); % ���ƹ�ת�ǣ�roll��  
legend('ekf', 'acc'); 
title('Roll Angle (degrees)');  
xlabel('Time');  
ylabel('Angle'); 


% ״̬ת�ƾ������
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

% �����������
function H_Update = update_H(q)
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    H_Update = [-2*q2  2*q3 -2*q0 2*q1;
                 2*q1  2*q0  2*q3 2*q2;
                 2*q0 -2*q1 -2*q2 2*q3];
end



