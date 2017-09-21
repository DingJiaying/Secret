%%测试M对估计性能的影响
%%整理之后的CRPHD针对自相关序列
clc;
close all;
clear all;

%% 信号参数
N = 200;
n = 1:1:N;
L = 100;

A =1;
B = 2;
phiA = 0.1*pi;
phiB = 0.2*pi;
w0 = 0.2*pi;   %从pi/25到pi-pi/25每隔pi/25取一个

sn = A*exp(- 1j*phiA)*exp(1j*w0*n) +  B*exp(- 1j*phiB)*exp(-1j*w0*n);

LL = 10:10:N-1;

%% 独立运行
T = 1000;
SNR = 20;
index = 1;
for L = LL
    wAcu = w0 * ones(1,T);             %真实的信号
    for t = 1:T
        xn = awgn(sn,SNR);
        for k=1:N-1
            rn(k)=r(xn,N,k);
        end
        Rn = rn - conj(rn);  %%新构造的序列
        LRn = length(Rn);  %%新构造的序列的长度
    %% SCRPHD方法的参数    
     Fre_SCRPHD = zeros(LRn,1); %%存储值
     for kk = LRn:-1:(L+2)      
         if kk > (L + 1)            
            SCRPHD_A = real(AR(Rn(kk:-1:(kk-L+1)), L));
            SCRPHD_B = real(BR(Rn(kk:-1:(kk-L+1)), L));
            argument_r = (SCRPHD_B+sqrt(SCRPHD_B^2+8*SCRPHD_A^2))/(4*SCRPHD_A);  
            if (argument_r>1)
                argument_r = 1;
            end
            if (argument_r<-1)
               argument_r = -1;
            end
         end
         Fre_SCRPHD(kk) = acos(argument_r);
     end
     SCRPHD_gu(t) = mean(Fre_SCRPHD(L+2:LRn));
         
    end
    MSE_SCRPHD(index) = 10*log10(mse(wAcu, SCRPHD_gu));
    index
    index = index + 1;
end
plot(LL, MSE_SCRPHD, '-^', 'LineWidth', 2);
hold on;

%仿真格式部分
legend('SCRPHD');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}L', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;