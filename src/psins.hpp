#pragma once

#define FASTLIO_PSINS_H
#ifdef FASTLIO_PSINS_H

#include <iostream>
#include <Eigen/Core>
#include "INSTools.h"

/**
 * @brief 
 * 利用前2周期的角速度、比力进行圆锥误差、划桨效应补偿
 * 要求IMU数据等间隔
 * 
 */
class psins
{
private:

    /* data */
public:
    bool                    cone_initd;     //圆锥误差初始化标志
    Eigen::Vector3d         cone;           //圆锥误差补偿量
    Eigen::Vector3d         scull;          //划桨效应补偿量
    Eigen::Vector3d         rot;            //旋转效应补偿量
    Eigen::Vector3d         Vcor;           //有害加速度补偿量
    
    double                  dt;             //时间间隔
    Eigen::Vector3d         wib[4];         //角速度
    Eigen::Vector3d         fib[4];         //比力
    Eigen::Vector3d         wt[3];          //角速度 * 时间
    Eigen::Vector3d         ft[3];          //比力 * 时间

    Eigen::Vector3d         wib_bia;        //角速度零偏
    Eigen::Vector3d         fib_bia;        //比力零偏

    Eigen::Vector3d         theta[3];       //姿态增量(b系)

    Eigen::Vector3d         delta_V_b;      //速度增量(b系)
    Eigen::Vector3d         delta_V_sf_n;   //旋转划桨效应补偿量 n系

    Eigen::Matrix3d         Cbn;            //当前姿态
    Eigen::Matrix3d         Cbn_pre;        //预估姿态

    Eigen::Vector3d         v_eb_n;         //当前速度
    Eigen::Vector3d         v_eb_n_pre;     //预估速度

    Eigen::Vector3d         p_eb_n;         //当前位置
    Eigen::Vector3d         p_eb_n_pre;     //预估位置
     

    psins();
    ~psins();

    bool INS_Data_Ready();
    bool INS_Cone();
    bool INS_Att_Predict();
    bool INS_RotScull();
    bool INS_Pos_Predict();
    bool INS_Start();
};

psins::psins(/* args */)
{
    dt = 0.005;
    cone_initd = false;

    wib_bia << 0,0,0;
    fib_bia << 0,0,0;
}

psins::~psins()
{
    cone_initd = true;
}

/**
 * @brief ins解算数据准备
 * 
 */
bool psins::INS_Data_Ready()
{
   int i;

    for(i=0; i<4; i++) 
    {
        wib[i] = wib[i] - wib_bia;
        fib[i] = fib[i] - fib_bia;
    }

    for(i=0; i<3; i++)
    {
         wt[i] = 0.5 * (wib[i] + wib[i+1]) * dt;
         ft[i] = 0.5 * (fib[i] + fib[i+1]) * dt;
    }

    return true;
}

/**
 * @brief 重叠三子样圆锥误差补偿
 * 
 */
bool psins::INS_Cone()
{
    double k_2_1 = 911.0 / 394;
    double k_1_0 = 911.0 / 394;
    double k_2_0 = 105.0 / 4998;
    double w_P = 89.0 / 54409;
    double w_Q = 14.0 / 249153;

    if(!cone_initd)
    {
        theta[2] = wt[2];
        theta[1] = wt[1];
        theta[0] = wt[0];
        cone_initd = true;
    }
    else
    {
        theta[2] = theta[1];
        theta[1] = theta[0];
        theta[0] = wt[0];
    }

    cone =  k_2_1 * wt[2].cross(wt[1]) +        \
            k_2_0 * wt[2].cross(wt[0]) +        \
            k_1_0 * wt[1].cross(wt[0]) +        \
            w_P * theta[1].cross(theta[0]) +    \
            w_Q * theta[2].cross(theta[0]);
    

    theta[0] = theta[0] + cone;

    return true;
}


/**
 * @brief 姿态积分
 * 
 */
bool psins::INS_Att_Predict()
{
    Eigen::Matrix3d     C_bp_bs;
    Eigen::Matrix3d     I3 = Eigen::Matrix3d::Identity();

    double mag_alpha = theta[0].norm();
    Eigen::Matrix3d  Alpha_ib_b = Askew(theta[0]); 

    double sin_alpha = sin(mag_alpha);
    double cos_alpha = sin(mag_alpha);
    double mag_alpha2 = mag_alpha * mag_alpha;

    if(mag_alpha > 1.0e-8)
    {
        double tmp1 = sin_alpha / mag_alpha;
        double tmp2 = (1.0 - cos_alpha) / mag_alpha;
        C_bp_bs = I3 + tmp1 * Alpha_ib_b + tmp2 * Alpha_ib_b * Alpha_ib_b;
    }
    else
    {
        C_bp_bs = I3 + Alpha_ib_b;
    }

    Cbn_pre = Cbn * C_bp_bs;

    return true;
}

/**
 * @brief 旋转划桨效应补偿
 * 
 */
bool psins::INS_RotScull()
{
    //旋转效应
    double alpha = wt[0].norm();
    double sin_alpha = sin(alpha);
    double cos_alpha = cos(alpha);
    double alpha2 = alpha * alpha;

    Eigen::Vector3d tmpwft = wt[0].cross(ft[0]);
    Eigen::Vector3d tmp1 = (1.0 - cos_alpha) / alpha2 * tmpwft;
    Eigen::Vector3d tmp2 = 1.0 / alpha2 * (1.0 - sin_alpha/alpha) * wt[0].cross(tmpwft);
    rot = tmp1 + tmp2;

    //划桨效应
    scull = 1.0 / 12 * (wt[1].cross(ft[0]) + ft[1].cross(wt[0]));

    //角速度曲线拟合
    Eigen::Vector3d a, b, c, dw, da, daw;
    Eigen::Vector3d w1, w2, w3, f1, f2, f3;
    double T = dt * 2;
    w1 = wib[2]; w2 = wib[1]; w3 = wib[0];
    a = w1;
    b = (4*w2 - w3 - 3*w1) / (2*T*T);
    c = 2 * (w1 - 2*w2 + w3) / (3*T*T);
    dw = 2*b + 6*c*T/2;

    //比力曲线拟合
    f1 = fib[2]; f2 = fib[1]; f3 = fib[0];
    a = f1;
    b = (4*f2 - f3 - 3*f1) / (2*T*T);
    c = 2 * (f1 - 2*f2 + f3) / (3*T*T);
    da = 2*b + 6*c*T/2;

    //剩余误差
    daw = 1/12 * da.cross(wib[1] + dw.cross(fib[1])) * dt * dt * dt;

    delta_V_b  = ft[0] + rot + scull - daw;

    return true;
}

/**
 * @brief 位置积分
 * 
 */
bool psins::INS_Pos_Predict()
{
    delta_V_sf_n = Cbn * delta_V_b;
    Eigen::Vector3d deltaV;
    deltaV = delta_V_sf_n + Vcor;

    p_eb_n_pre = p_eb_n + 0.5 * (v_eb_n + v_eb_n_pre);
    return true;
}

/**
 * @brief 位置积分
 * 
 */
bool psins::INS_Start()
{
    //圆锥误差补偿
    INS_Cone();

    //姿态积分
    INS_Att_Predict();

    //旋转划桨效应补偿
    INS_RotScull();

    //位置积分
    INS_Pos_Predict();

    return true;
}

#endif