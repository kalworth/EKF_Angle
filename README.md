前言
    浅浅记录一下对KF和四元数姿态解算的学习，欢迎大家指正。对于EKF和四元数，网上的资料已经很充分了，DR_CAN的公式推导真的很详细，这里不再赘述，该文对 
𝑤
𝑏
𝑖
𝑎
𝑠
 没有更新，理解原理之后自行加上即可。
基于EKF的姿态解算 - 知乎 (zhihu.com)
四元数EKF姿态更新算法 - 知乎 (zhihu.com)
使用扩展卡尔曼滤波（EKF）进行AHRS九轴姿态融合_ahrs 卡尔曼滤波_@奔跑的蜗牛@的博客-CSDN博客
【卡尔曼滤波器】1_递归算法_Recursive Processing_哔哩哔哩_bilibili
实现步骤
1.状态量(先验估计)
    在KF中，这里设置的状态量 
𝑥
𝑘
 即是机体的四元数 
(4)
𝑞
=
[
𝑞
0
𝑞
1
𝑞
2
𝑞
3
]
𝑇
 
(5)
𝑥
𝑘
=
𝑞
 
    状态方程如下
(6)
𝑥
𝑘
=
𝐴
𝑥
𝑘
−
1
+
𝑊
     其中 
𝑊
 ~ 
𝑁
(
0
,
𝑄
)
 ,我们实际能够计算的只有前半部分
(7)
𝑥
^
𝑘
−
=
𝐴
𝑥
^
𝑘
−
1
 
    状态转移矩阵A要从四元数微分方程中得来
    已知四元数微分方程如下:
(8)
𝑞
˙
=
1
2
[
0
−
𝑤
𝑥
−
𝑤
𝑦
−
𝑤
𝑧
𝑤
𝑥
0
𝑤
𝑧
−
𝑤
𝑦
𝑤
𝑦
−
𝑤
𝑧
0
𝑤
𝑥
𝑤
𝑧
𝑤
𝑦
−
𝑤
𝑥
0
]
𝑞
 
(9)
Ω
=
[
0
−
𝑤
𝑥
−
𝑤
𝑦
−
𝑤
𝑧
𝑤
𝑥
0
𝑤
𝑧
−
𝑤
𝑦
𝑤
𝑦
−
𝑤
𝑧
0
𝑤
𝑥
𝑤
𝑧
𝑤
𝑦
−
𝑤
𝑥
0
]
 
    在数字系统中，更新四元数的方式
(10)
𝑞
(
𝑡
+
Δ
𝑡
)
=
(
𝐼
4
+
1
2
Ω
Δ
𝑡
)
𝑞
𝑡
     即得完整的状态方程
(11)
𝑥
^
𝑘
−
=
(
𝐼
4
+
1
2
Ω
Δ
𝑡
)
𝑥
^
𝑘
−
1
     这里，我们得到了状态转移矩阵A
(12)
𝐴
=
(
𝐼
4
+
1
2
Ω
Δ
𝑡
)
 
2.计算先验误差协方差
(13)
𝑃
𝑘
−
=
𝐴
𝑃
𝑘
−
1
𝐴
𝑇
+
𝑄
     这里套公式即可。
3.观测量
    观测方程如下
(14)
𝑍
𝑘
=
𝐻
𝑥
𝑘
+
𝑉
 
其中 
𝑉
 ~ 
𝑁
(
0
,
𝑅
)
 
这里观测量 
𝑍
𝑘
 就是加速度计采集到的数据。
(15)
𝑍
𝑘
=
[
𝑎
𝑥
𝑎
𝑦
𝑎
𝑧
]
𝑇
 在东北天坐标系下，我们定义重力加速度 
𝐺
 向量，以及机体的加速度向量 
𝐺
𝑏
 
(16)
𝐺
=
[
0
0
1
]
𝑇
 
(17)
𝐺
𝑏
=
[
𝑎
𝑥
′
𝑎
𝑦
′
𝑎
𝑧
′
]
𝑇
 根据 
𝐺
=
𝐶
𝑛
𝑏
𝐺
𝑏
 
(18)
𝐺
𝑏
=
𝐶
𝑏
𝑛
[
0
0
1
]
=
[
2
𝑞
1
𝑞
3
−
2
𝑞
0
𝑞
2
2
𝑞
2
𝑞
3
+
2
𝑞
0
𝑞
1
1
−
2
𝑞
1
2
−
2
𝑞
2
2
]
 前面的文章中里 
1
−
2
𝑞
1
2
−
2
𝑞
2
2
 与 
𝑞
0
2
−
𝑞
1
2
−
𝑞
2
2
+
𝑞
3
2
 等价，这里使用的四元数都是单位四元数，在具体代码中一定要归一化。对上式求雅克比矩阵，我们可以得到观测矩阵 
𝐻
 
(19)
𝐻
=
[
−
2
𝑞
2
2
𝑞
3
−
2
𝑞
0
2
𝑞
1
2
𝑞
1
2
𝑞
0
2
𝑞
3
2
𝑞
2
2
𝑞
0
−
2
𝑞
1
−
2
𝑞
2
2
𝑞
3
]
 4.计算卡尔曼增益
(20)
𝐾
𝑘
=
𝑃
𝑘
−
𝐻
𝑇
𝐻
𝑃
𝑘
−
𝐻
𝑇
+
𝑅
 这里也套公式。
5.后验估计
(21)
𝑥
^
𝑘
=
𝑥
^
𝑘
−
1
−
+
𝐾
𝑘
(
𝑍
𝑘
−
𝐻
𝑥
^
𝑘
−
1
−
)
 6.更新误差协方差
(22)
𝑃
𝑘
=
(
𝐼
4
−
𝐾
𝑘
𝐻
)
𝑃
𝑘
−
 整个计算流程就是这些了。在得到观测矩阵 
𝐻
 中，把原本的非线性函数(18)，通过求雅可比矩阵，相当于一级的泰勒展开，完成了线性化。
Python源码
github地址：kalworth/EKF.py: 基于扩展卡尔曼(EKF)的四元数姿态解算(六轴) (github.com)
使用的数据是通过串口打印保存下来的，顺序是ax,ay,az,gx,gy,gz
2024/4/29注：python代码运行的时候，注意做单位转换，xlxs里面的角速度单位是deg/s，加速度单位是g，需要转换成rad/s以及m²/s再进行运算，注意！！！(当前python代码里没有加，懒得改了，后面有C和Matlab版本的)

加速度计与陀螺仪数据
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

imu_data = pd.read_csv("C:/Serial Debug 2023-8-12 174233.csv")
imu_data = np.asarray(imu_data)

Q = [0.01] * 4
R = [1000.0] * 3
I = [1.0] * 4
P = [10000.0] * 4

class EKF():
    def __init__(self, period):
        self.Q_matrix  = np.diag(Q)
        self.R_matrix  = np.diag(R)
        self.I_matrix  = np.diag(I)
        self.halfT     = 1/2 * period
        self.P_matrix  = np.diag(P)
        self.A_matrix  = np.zeros([4,4])
        self.H_matrix  = np.zeros([3,4])
        self.K_matrix  = np.zeros([4,3])
        self.K_vector  = np.zeros([4,])
        self.T_vector  = np.zeros([3,])
        self.HX_vector = np.zeros([3,])
        self.Z_vector  = np.zeros([3,])
        self.pitch     = 0.0
        self.roll      = 0.0
        #self.q         = np.random.randn(4)
        self.q         = np.array([1.0,0.0,0.0,0.0])
        self.pitch_list= []
        self.roll_list = []
        self.a_pitch   = 0.0
        self.a_roll    = 0.0
        self.a_pitch_list = []
        self.a_roll_list = []

    def normalizeQuternion(self, q: np.ndarray):

        norm = np.linalg.norm(q,2)

        norm_q = q/norm

        return  norm_q

    def priori(self, gx, gy, gz):

        gx_ = gx * self.halfT
        gy_ = gy * self.halfT
        gz_ = gz * self.halfT

        self.A_matrix = np.array([
            [1, -gx_, -gy_, -gz_],
            [gx_, 1,   gz_, -gy_],
            [gy_, -gz_, 1,   gx_],
            [gz_, gy_, -gx_,   1]
        ])

        self.q = np.dot(self.A_matrix,self.q)

    def cal_P_matrix(self):
        self.P_matrix = np.dot(np.dot(self.A_matrix,self.P_matrix),np.transpose(self.A_matrix)) + self.Q_matrix

    def cal_HX_vector(self, q,ax, ay, az):
        self.H_matrix = np.array([
            [-2*q[2], 2*q[3], -2*q[0], 2*q[1]],
            [2*q[1],  2*q[0],  2*q[3], 2*q[2]],
            [2*q[0], 2*q[1], 2*q[2], 2*q[3]]
        ])
        self.HX_vector = np.array([
            2*q[1]*q[3] - 2*q[0]*q[2],
            2*q[2]*q[3] + 2*q[0]*q[1],
            1 - 2*(q[1]**2) - 2*(q[2]**2)
        ])

        self.Z_vector = np.array([ax,ay,az])

    def cal_K_matrix(self):
        self.K_matrix = np.dot(
            np.dot(self.P_matrix,np.transpose(self.H_matrix)),
            np.linalg.inv(
                np.dot(
                    np.dot(
                        self.H_matrix,self.P_matrix
                    ),np.transpose(self.H_matrix)
                ) + self.R_matrix
            )
        )

    def posterior(self):
        # 原来是self.q_k 改正为 self.q
        self.q = self.q + np.dot(
            self.K_matrix,
            self.Z_vector - self.HX_vector
        )

    def update_P_matrix(self):
        self.P_matrix = np.dot(
            self.I_matrix - np.dot(self.K_matrix,self.H_matrix),
            self.P_matrix
        )

    def Quternion2Angle(self,q):
        self.roll = -np.arcsin(2*(q[0]*q[2] - q[1]*q[3])) * 57.3
        self.pitch = np.arctan2((2*(q[0]*q[1] + q[1]*q[3])), (2*(q[0]*q[0] + q[3]*q[3]) - 1.0)) * 57.3
        self.pitch_list.append(self.pitch)
        self.roll_list.append(self.roll)

    def plot_angle(self, angle_list, color, angle_name):
         x = np.arange(0, len(angle_list), 1)
         y = np.array(angle_list)
         plt.plot(x,y,color=color,linewidth=1,label=angle_name)

    def cal_a_angle(self, ax, ay, az):
        self.a_pitch = np.arctan2(ay,az) * 57.3
        self.a_roll = np.arctan(ax/np.sqrt(ay*ay+az*az)) * 57.3
        self.a_pitch_list.append(self.a_pitch)
        self.a_roll_list.append(self.a_roll)


    def EKF_update(self,imu_sorce_data):
        for i in range(2):
            for data in imu_sorce_data:
                self.q = self.normalizeQuternion(self.q)
                self.priori(data[3], data[4], data[5])
                self.q = self.normalizeQuternion(self.q)
                self.cal_P_matrix()
                self.cal_HX_vector(self.q, data[0], data[1], data[2])
                self.cal_K_matrix()
                self.posterior()
                self.update_P_matrix()
                self.Quternion2Angle(self.q)
                self.cal_a_angle(data[0],data[1],data[2])

        plt.figure()  # 画布尺寸默认
        self.plot_angle(self.pitch_list,'red','pitch')
        self.plot_angle(self.a_pitch_list, 'blue', 'roll')
        plt.legend(['EKF_pitch','acc_pitch'], loc='best')
        plt.figure()  # 画布尺寸默认
        self.plot_angle(self.roll_list,'green','pitch')
        self.plot_angle(self.a_roll_list, 'black', 'roll')
        plt.legend(['EKF_roll','acc_roll'], loc='best')
        plt.show()

if __name__ == "__main__":
    EKF_model = EKF(period=0.0001)
    EKF_model.EKF_update(imu_data)
   对滤波效果影响较大的因素：机体四元数q的初值以及采样周期。卡尔曼增益会收敛，QR阵参数合理即可。

roll角经过EKF与只用加速度计解算作比较

局部放大
    比较明显的是经过EKF之后，波形更加平滑，利于进一步的使用。
    以下是在MCU中运行的程序，如果有错的地方欢迎指出来。
typedef struct {
    float period;                       // 姿态更新周期
    float half_T;                       // 姿态更新周期的1/2      单位：s

    float I_matrix[16];                 // 单位矩阵              4 * 4
    float Q_matrix[16];                 // 过程噪声协方差矩阵      4 * 4
    float R_matrix[9];                  // 观测噪声协方差矩阵      3 * 3
    float A[16];                        // 状态转移矩阵           4 * 4
    float H[12];                        // 观测矩阵              3 * 4
    float p_hat[16];                    // 误差协方差矩阵         4 * 4
    float K[12];                        // 卡尔曼增益矩阵         4 * 3
    float Kq[4];                        // 四元数卡尔曼增益        4 * 1
    float T[3];                         // 残差向量               3 * 1
    float q[4];                         // 机体四元数             4 * 1
    float q_k[4];                       // 优化后的四元数          4 * 1

    float vector_hat[3];                // 观测重力加速度          3 * 1

    void (*filter)(float,float,float,   // EKF更新姿态函数指针
                        float,float,float);
}EKF_StructTypeDef;

    关于EKF的结构体
static EKF_StructTypeDef      imu_ekf;                              // EKF
     然后是初始化，  噪声矩阵与我们的状态量和观测量是对应的。 协方差矩阵对应预测向量的协方差，  协方差矩阵对应观测向量的协方差。  的值与预测值的权重成反比，即Q矩阵中协方差越小，越相信预测值(对应相信角速度)，R的值与观测值的权重成反比，即R矩阵中协方差越小，越相信观测值(对应相信加速度)。
/**
  * @brief  EKF初始化
  * @param  None
  * @note   None
  * @retval None
  */
static void EKF_init(void)
{

    imu_ekf.period = 0.001;  // 姿态解算周期 1ms
    imu_ekf.half_T = imu_ekf.period / 2.0f;

    imu_ekf.I_matrix[0] = 1; imu_ekf.I_matrix[1] = 0;   imu_ekf.I_matrix[2] = 0;  imu_ekf.I_matrix[3] = 0;
    imu_ekf.I_matrix[4] = 0; imu_ekf.I_matrix[5] = 1;   imu_ekf.I_matrix[6] = 0;  imu_ekf.I_matrix[7] = 0;
    imu_ekf.I_matrix[8] = 0; imu_ekf.I_matrix[9] = 0;   imu_ekf.I_matrix[10] = 1; imu_ekf.I_matrix[11] = 0;
    imu_ekf.I_matrix[12] = 0;imu_ekf.I_matrix[13] = 0;  imu_ekf.I_matrix[14] = 0; imu_ekf.I_matrix[15] = 1;

    float Q_Val = 0.01;
    imu_ekf.Q_matrix[0] = Q_Val; imu_ekf.Q_matrix[1] = 0;     imu_ekf.Q_matrix[2] = 0;     imu_ekf.Q_matrix[3] = 0;
    imu_ekf.Q_matrix[4] = 0;     imu_ekf.Q_matrix[5] = Q_Val; imu_ekf.Q_matrix[6] = 0;     imu_ekf.Q_matrix[7] = 0;
    imu_ekf.Q_matrix[8] = 0;     imu_ekf.Q_matrix[9] = 0;     imu_ekf.Q_matrix[10] = Q_Val;imu_ekf.Q_matrix[11] = 0;
    imu_ekf.Q_matrix[12] = 0;    imu_ekf.Q_matrix[13] = 0;    imu_ekf.Q_matrix[14] = 0;    imu_ekf.Q_matrix[15] = Q_Val;

    float R_Val = 1000000;
    imu_ekf.R_matrix[0] = R_Val; imu_ekf.R_matrix[1] = 0;     imu_ekf.R_matrix[2] = 0;
    imu_ekf.R_matrix[3] = 0;     imu_ekf.R_matrix[4] = R_Val; imu_ekf.R_matrix[5] = 0;
    imu_ekf.R_matrix[6] = 0;     imu_ekf.R_matrix[7] = 0;     imu_ekf.R_matrix[8] = R_Val;

    float P_Val = 100000;
    imu_ekf.p_hat[0] = P_Val;imu_ekf.p_hat[1] = 0;    imu_ekf.p_hat[2] = 0;     imu_ekf.p_hat[3] = 0;
    imu_ekf.p_hat[4] = 0;    imu_ekf.p_hat[5] = P_Val;imu_ekf.p_hat[6] = 0;     imu_ekf.p_hat[7] = 0;
    imu_ekf.p_hat[8] = 0;    imu_ekf.p_hat[9] = 0;    imu_ekf.p_hat[10] = P_Val;imu_ekf.p_hat[11] = 0;
    imu_ekf.p_hat[12] = 0;   imu_ekf.p_hat[13] = 0;   imu_ekf.p_hat[14] = 0;    imu_ekf.p_hat[15] = P_Val;

    imu_ekf.q[0] = 1;
    imu_ekf.q[1] = 0;
    imu_ekf.q[2] = 0;
    imu_ekf.q[3] = 0;

    imu_ekf.q_k[0] = 1;
    imu_ekf.q_k[1] = 0;
    imu_ekf.q_k[2] = 0;
    imu_ekf.q_k[3] = 0;

    imu_ekf.filter = ekf_update;
}
接下来就是四元数的更新了，实现上文提到的过的完整的EKF的计算流程，关于矩阵运算的函数，我找的GPT生成的，自己用dev-c还验证了一下，计算卡尔曼增益的时候需要求矩阵的逆阵，我只用了一个只能求3×3的逆矩阵的函数。(事实是，没找到合适的矩阵库，参考的几位大佬用的是一维的数组，设定一个形状的矩阵)。想要获取欧拉角，做一步转换即可。
/**
  * @brief  扩展卡尔曼 优化四元数
  * @param  gx      x轴角速度
  * @param  gy      y轴角速度
  * @param  gz      z轴角速度
  * @param  ax      x轴加速度
  * @param  ay      y轴加速度
  * @param  az      z轴加速度
  * @note   None
  * @retval None
  */
static void ekf_update(float gx, float gy, float gz, float ax, float ay, float az)
{
#define DIV_180_PI       (float)(57.295779513f)
#define GRAVITY          (float)(9.8f)

    // 单位换算, deg/s转换到rad/s, g换算到m²/2
    gx /= DIV_180_PI;
    gy /= DIV_180_PI;
    gz /= DIV_180_PI;

    ax *= GRAVITY;
    ay *= GRAVITY;
    az *= GRAVITY;

    normalize_quternion(imu_ekf.q);
    /****************************** 状态方程 ******************************/
    /* 1/2姿态更新周期赋值 */
    float half_T = imu_ekf.half_T;
    /* 四元数微分方程离散化 */
    imu_ekf.q[0] += (-imu_ekf.q[1] * gx - imu_ekf.q[2] * gy - imu_ekf.q[3] * gz) * (float)half_T;
    imu_ekf.q[1] += ( imu_ekf.q[0] * gx + imu_ekf.q[2] * gz - imu_ekf.q[3] * gy) * (float)half_T;
    imu_ekf.q[2] += ( imu_ekf.q[0] * gy - imu_ekf.q[1] * gz + imu_ekf.q[3] * gx) * (float)half_T;
    imu_ekf.q[3] += ( imu_ekf.q[0] * gz + imu_ekf.q[1] * gy - imu_ekf.q[2] * gx) * (float)half_T;
    /*   四元数 归一化    */
    normalize_quternion(imu_ekf.q);
    /*  状态转移矩阵赋值   */
    imu_ekf.A[0] =            1;  imu_ekf.A[1] = -gx * half_T; imu_ekf.A[2] = -gy * half_T;  imu_ekf.A[3] = -gz * half_T;
    imu_ekf.A[4] =  gx * half_T;  imu_ekf.A[5] =            1; imu_ekf.A[6] =  gz * half_T;  imu_ekf.A[7] = -gy * half_T;
    imu_ekf.A[8] =  gy * half_T;  imu_ekf.A[9] = -gz * half_T; imu_ekf.A[10]=            1;  imu_ekf.A[11]=  gx * half_T;
    imu_ekf.A[12]=  gz * half_T;  imu_ekf.A[13]=  gy * half_T; imu_ekf.A[14]= -gx * half_T;  imu_ekf.A[15]=            1;
    /****************************** 状态方程 ******************************/

    /**************************** 估算协方差矩阵 ***************************/
    float temp_matrix1[16];
    float temp_matrix2[16];
    float temp_matrix3[16];
    matrix_multiply(imu_ekf.A,imu_ekf.p_hat,temp_matrix1,4,4,4);
    matrix_transpose(imu_ekf.A,temp_matrix2,4,4);
    matrix_multiply(temp_matrix1,temp_matrix2,temp_matrix3,4,4,4);
    matrix_add(temp_matrix3,imu_ekf.Q_matrix,imu_ekf.p_hat,4,4);
    /**************************** 估算协方差矩阵 ***************************/

    /****************************** 观测方程 ******************************/
    /*   观测向量计算  */
    imu_ekf.vector_hat[0] = 2 * (imu_ekf.q[1] * imu_ekf.q[3] - imu_ekf.q[0] * imu_ekf.q[2]);
    imu_ekf.vector_hat[1] = 2 * (imu_ekf.q[2] * imu_ekf.q[3] + imu_ekf.q[0] * imu_ekf.q[1]);
    imu_ekf.vector_hat[2] = imu_ekf.q[0] * imu_ekf.q[0] - imu_ekf.q[1] * imu_ekf.q[1] - imu_ekf.q[2] * imu_ekf.q[2] + imu_ekf.q[3] * imu_ekf.q[3];
    /*    观测矩阵赋值    */
    imu_ekf.H[0] = -2 * imu_ekf.q[2]; imu_ekf.H[1] =  2 * imu_ekf.q[3]; imu_ekf.H[2] =  -2 * imu_ekf.q[0]; imu_ekf.H[3] =  2 * imu_ekf.q[1];
    imu_ekf.H[4] =  2 * imu_ekf.q[1]; imu_ekf.H[5] =  2 * imu_ekf.q[0]; imu_ekf.H[6] =   2 * imu_ekf.q[3]; imu_ekf.H[7] =  2 * imu_ekf.q[2];
    imu_ekf.H[8] =  2 * imu_ekf.q[0]; imu_ekf.H[9] = -2 * imu_ekf.q[1]; imu_ekf.H[10] = -2 * imu_ekf.q[2]; imu_ekf.H[11] = 2 * imu_ekf.q[3];
    /****************************** 观测方程 ******************************/

    /**************************** 计算卡尔曼增益 ***************************/
    float temp_matrix4[12];
    float temp_matrix5[12];
    float temp_matrix6[9];
    float temp_matrix7[9];
    float temp_matrix8[9];
    float HT[12];
    matrix_transpose(imu_ekf.H,HT,3,4);
    matrix_multiply(imu_ekf.p_hat,HT,temp_matrix4,4,4,3);
    matrix_multiply(imu_ekf.H,imu_ekf.p_hat,temp_matrix5,3,4,4);
    matrix_multiply(temp_matrix5,HT,temp_matrix6,3,4,3);
    matrix_add(temp_matrix6,imu_ekf.R_matrix,temp_matrix7,3,3);
    calculateInverse(temp_matrix7,temp_matrix8);
    matrix_multiply(temp_matrix4,temp_matrix8,imu_ekf.K,4,3,3);
    /**************************** 计算卡尔曼增益 ***************************/

    /****************************** 计算残差 ******************************/
    /*   归一化加速度计    */
    float norm = sqrtf(ax * ax + ay * ay + az* az);
    ax /= norm;
    ay /= norm;
    az /= norm;
    /*    计算残差向量    */
    imu_ekf.T[0] = ax - imu_ekf.vector_hat[0];
    imu_ekf.T[1] = ay - imu_ekf.vector_hat[1];
    imu_ekf.T[2] = az - imu_ekf.vector_hat[2];
    /****************************** 计算残差 ******************************/

    /*************************** 后验估计四元数 ****************************/
    matrix_multiply(imu_ekf.K,imu_ekf.T,imu_ekf.Kq,4,3,1);
    /*    后验四元数      */
    imu_ekf.q_k[0] = imu_ekf.q[0] + imu_ekf.Kq[0];
    imu_ekf.q_k[1] = imu_ekf.q[1] + imu_ekf.Kq[1];
    imu_ekf.q_k[2] = imu_ekf.q[2] + imu_ekf.Kq[2];
    imu_ekf.q_k[3] = imu_ekf.q[3] + imu_ekf.Kq[3];
    /*    归一化四元数     */
    normalize_quternion(imu_ekf.q_k);
    /*    更新 四元数     */
    imu_ekf.q[0] = imu_ekf.q_k[0];
    imu_ekf.q[1] = imu_ekf.q_k[1];
    imu_ekf.q[2] = imu_ekf.q_k[2];
    imu_ekf.q[3] = imu_ekf.q_k[3];
    /*************************** 后验估计四元数 ****************************/

    /*************************** 更新协方差矩阵 ****************************/
    float temp_matrix9[16];
    float temp_matrix10[16];
    float temp_matrix11[16];
    matrix_multiply(imu_ekf.K,imu_ekf.H,temp_matrix9,4,3,4);
    matrix_subtract(imu_ekf.I_matrix,temp_matrix9,temp_matrix10,4,4);
    matrix_multiply(temp_matrix10,imu_ekf.p_hat,temp_matrix11,4,4,4);
    matrix_copy(temp_matrix11,imu_ekf.p_hat,4,4);
    /*************************** 更新协方差矩阵 ****************************/
}
重要的事情说三遍，一定要注意采样周期！一定要注意采样周期！一定要注意采样周期！

最后，结尾留一个问题，这里的q的初始值为什么是  ,怪中之怪。(见评论区)
优化效果
这里和纯加速度做对比。

加速度计解算与ekf解算

局部放大
Vofa+上位机查看mcu滤波效果

红色曲线为EKF优化之后的俯仰角，蓝色曲线为纯加速度计解算
这幅图是大角度小摆动速度和小角度大摆动速度测得的解算曲线，很明显的可以看到，EKF对运动加速度带来的误差有不错的抑制效果。
姿态解算修正
今天在做matlab仿真时，发现姿态更新的单位有点问题，上面的C代码部分已经同步更新。

Roll角

Pitch角
在修改了单位之后，可以对比之前的图片，相位滞后几乎没有了，再用上玺佬的PQR三个矩阵的初始值，效果可以在平衡车上正常使用了。
最近比较忙，矩阵运算部分过一阵子会开源出来，不限制MCU的平台，纯软件的矩阵算法，因为这个EKF只用到了特殊尺寸的矩阵运算，所以库里有针对性的矩阵运算。
MATLAB仿真程序以及串口输出文件
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
eul_acc_deg = zeros(numRows*epoch,2);

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
R_value = 1000000;
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
        q = q + K * (Z - H*q);

        % 更新误差协方差
        P = (eye(4) - K * H) * P;

        % 欧拉角解算
        qt = transpose(q);
        eul = quat2eul(qt);  
        eul_deg(deg_num,:) = rad2deg(eul);

        eul_acc_deg(deg_num,1) = rad2deg(atan2(ay,az));
        eul_acc_deg(deg_num,2) = -rad2deg(atan(ax/sqrt(ay*ay+az*az)));
        deg_num = deg_num + 1;
    end
end

line_width = 2;
figure; % 激活第三个子图  
plot(1:numRows*epoch, eul_deg(:,3), 'b-','LineWidth', line_width); % 绘制滚转角（roll）  
hold on; % 保持当前图形，以便添加更多曲线  
plot(1:numRows*epoch, eul_acc_deg(:,1), 'r-','LineWidth', line_width); % 绘制滚转角（roll）  
legend('ekf', 'acc');  
title('Pitch Angle (degrees)');  
xlabel('Time');  
ylabel('Angle'); 

figure; % 激活第三个子图  
plot(1:numRows*epoch, eul_deg(:,2), 'b-','LineWidth', line_width); % 绘制滚转角（roll）  
hold on; % 保持当前图形，以便添加更多曲线  
plot(1:numRows*epoch, eul_acc_deg(:,2), 'r-','LineWidth', line_width); % 绘制滚转角（roll）  
legend('ekf', 'acc'); 
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
matlab仿真脚本以及xlxs串口数据地址kalworth/EKF_Angle (github.com)
