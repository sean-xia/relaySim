# 中继仿真代码的说明

*本篇是关于A New Look at Dual-Hop Relaying: Performance Limits with Hardware Impairments论文仿真代码的说明*

```matlab
rng('shuffle'); 
```


设置产生随机数的种子与系统时间有关系，这样每次产生的随机数都不不可能一样。

```matlab
nbrOfIterations = 1000000;
```
设置迭代次数，为了和解析结果拟合的足够好，此处设置的信道产生次数应尽可能多。

```matlab
scenario = 2;
```
本例计算了两种情形，一种是瑞利衰落信道1，另一个是Nakagami信道2，此处假设为第二种情况，当然可以修改为瑞利信道，将其值设置为1

```matlab
if scenario == 1 %瑞利衰落信道情况
    alpha1 = 1; %对四种情况没有意义
    beta1 = 1;  %第一跳（S到R）的方差
    
    %第二跳信道的参数
    alpha2 = 1; %对于瑞利信道无意义
    beta2 = 1;  %第二跳的信道方差
    
   
    RHO = [beta1*exprnd(1,1,nbrOfIterations); beta2*exprnd(1,1,nbrOfIterations)];
 ```
此处$\rho$是随机信道增益的实现，产生了两跳信道的$h_i$, $\rho_i = |h_i|^2$. 由于瑞利衰落信道的增益的平方是指数分布，文章中假设信道增益的方差为1，那么信道增益的平方即是均值为1的指数分布。$\rho_i\sim exp(1)$

 ```matlab   
elseif scenario == 2 %Nakagami-m 衰落情形
  %第一跳的信道参数  
    alpha1 = 2; %alpha fading parameter at first hop
    beta1 = 1;  %beta fading parameter at first hop
  %第二跳的信道参数  
    alpha2 = 2; %alpha fading parameter at second hop
    beta2 = 1;  %beta fading parameter at second hop
    
    %产生衰落信道参数，这里rho是信道系数范数的平方
    RHO = [gamrnd(alpha1,beta1,1,nbrOfIterations); gamrnd(alpha2,beta2,1,nbrOfIterations)];
    
end
```
对于Nakagami信道而言，信道增益$\rho_i=|h_i|^2\sim \Gamma(\alpha_i,\beta_i)$ ,也就是服从Gamma分布。

```matlab
R = [2 5];   %传输的目标信息速率 (每两次信道使用)
x = 2.^R-1;  %对应的中断信噪比门限

N1 = 1;      %中继节点的规范化噪声方差
N2 = 1;      %目标节点的规范化噪声方差

SNR_dB = 0:5:40; %第一跳的信噪比范围
P1_dB = SNR_dB - 10*log10(alpha1*beta1); %相应的信源传输功率范围
P1 = 10.^(P1_dB/10); %转化成线性功率值

P2 = P1; %中继上采用同样的发射功率.
```
此处，信道功率增益为$E_{\rho_i}\{\rho_i\}=\alpha_i\beta_i$, 转化成分贝单位后和信噪比相减，得到发射功率。注意，信道增益功率为负数。


```matlab
kappa1 = 0.1; %第一跳收发机的误差向量幅度
kappa2 = 0.1; %第二跳收发机得误差向量幅度
d = kappa1^2+kappa2^2+kappa1^2*kappa2^2; %Recurring function of kappa1 and kappa2
```
$d$是论文中公式（13）和（14）中的表达式 $d=\kappa_1^2+\kappa2^2+\kappa1^2\kappa2^2$

```matlab

%蒙特卡洛仿真中断概率的结果存储变量初始化：放大转发情形
OP_ideal_af_f = zeros(length(x),length(P1)); %AF固定增益，理想硬件的中断概率
OP_ideal_af_v = zeros(length(x),length(P1)); %AF变增益，理想硬件的中断概率

OP_nonideal_af_f = zeros(length(x),length(P1)); %AF固定增益，非理想硬件的中断概率
OP_nonideal_af_v = zeros(length(x),length(P1)); %AF变增益，非理想硬件的中断概率


%蒙特卡洛仿真中断概率的结果存储变量初始化：译码转发情形
OP_ideal_df = zeros(length(x),length(P1));    %DF理想硬件的中断概率
OP_nonideal_df = zeros(length(x),length(P1)); %DF非理想硬件的中断概率

```

%用于计算中继的平均衰落功率$\mathbb{E}_{\rho_1}[\rho_1]=\alpha_1\beta_1$
```matlab
exp_rho1 = alpha1*beta1;
```

以下用蒙特卡洛仿真计算中断概率。假设
* 源节点发射功率变化以实现设定的平均信噪比
* 中继节点发射功率和源节点采用相同的变化方式

```matlab

for p = 1:length(P1)
    
    G_ideal_f = sqrt(P2(p)/(P1(p)*exp_rho1 + N1));     % (1)Ideal case Fixed Gain
    G_nonideal_f = sqrt(P2(p)/(P1(p)*exp_rho1*(1 + kappa1^2) + N1));   % (2)With impairments Fixed Gain
    

```
(1) 计算固定增益，理想收发机情况下的放大因子，即文中（15）式
$$G_{id}^f=\sqrt{\frac{P_2}{P_1\mathbb{E}_{\rho_1}\{\rho_1\}+N_1}}$$

(2)计算固定增益，非理想收发机情况下的放大因子，即文中（11）式
$$G_{ni}^f=\sqrt{\frac{P_2}{P_1\mathbb{E}_{\rho_1}\{\rho_1\}(1+\kappa_1^2)+N_1}}$$


```matlab
    %蒙特卡洛仿真
    for kk = 1:nbrOfIterations
        
        rho1 = RHO(1,kk);  %提取每次信道实现的信道增益
        rho2 = RHO(2,kk);
        
        %计算AF模式下的端到端信噪失真比（SNDR）
        SNDR_ideal_af_f = (rho1*rho2)/(rho2*N1/P1(p) + N2/G_ideal_f^2/P1(p)); %（3）Ideal hardware, fixed gain
        SNDR_ideal_af_v = (rho1*rho2)/(rho2*N1/P1(p) + N2*(rho1/P2(p)+N1/P2(p)/P1(p))); （4）%Ideal hardware, variable gain
        
        SNDR_nonideal_af_f = (rho1*rho2)/(rho1*rho2*d + N1*rho2/P1(p) *(1+kappa2^2) + N2/G_nonideal_f^2/P1(p)); %（5）Non-ideal hardware, fixed gain
        SNDR_nonideal_af_v = (rho1*rho2)/(rho1*rho2*d + rho1*N2*(1+kappa1^2)/P2(p) + rho2*N1*(1+kappa2^2)/P1(p) + N1*N2/P1(p)/P2(p));  %（6）Non-ideal hardware, variable gain
```
(3)计算理想硬件固定增益的SNDR，文中（16）式
(4)计算理想硬件变化增益的SNDR，文中（16）式
(5)计算非理想硬件，固定增益的SNDR，文中（13）式
(6)计算非理想硬件，变化增益的SNDR，文中（14）式
```matlab
        %Compute signal-to-noise-and-distortion ratios (SNDRs) with DF relaying
        SNDR_ideal_df = min([(rho1*P1(p))/N1,(rho2*P2(p))/N2]); %Ideal hardware  文中( 18)式
        SNDR_nonideal_df = min([(rho1*P1(p))/(P1(p)*rho1*kappa1^2 + N1),(rho2*P2(p))/(P2(p)*rho2*kappa2^2 + N2)]); %Non-ideal hardware文中（17）式
        
        
        %分别针对六种情况判断是否发生中断，（AF四种，DF两种）检测条件： (if SNDR < x). 
        %将中断次数存放到相应的变量中去。
        for ind = 1:length(x)
            
            if SNDR_ideal_af_f < x(ind) %Fixed gain AF relaying, ideal hardware
                OP_ideal_af_f(ind,p) = OP_ideal_af_f(ind,p) + 1;
            end
            
            if SNDR_nonideal_af_f < x(ind) %Fixed gain AF relaying, non-ideal hardware
                OP_nonideal_af_f(ind,p) = OP_nonideal_af_f(ind,p) + 1;
            end
            
            if SNDR_ideal_af_v < x(ind) %Variable gain AF relaying, ideal hardware
                OP_ideal_af_v(ind,p) = OP_ideal_af_v(ind,p) + 1;
            end
            
            if SNDR_nonideal_af_v < x(ind) %Variable gain AF relaying, non-ideal hardware
                OP_nonideal_af_v(ind,p) = OP_nonideal_af_v(ind,p) + 1;
            end
            
            if SNDR_ideal_df < x(ind) %DF relaying, ideal hardware
                OP_ideal_df(ind,p) = OP_ideal_df(ind,p) + 1;
            end
            
            if SNDR_nonideal_df < x(ind) %DF relaying, non-ideal hardware
                OP_nonideal_df(ind,p) = OP_nonideal_df(ind,p) + 1;
            end
        end
        
    end
end
```
