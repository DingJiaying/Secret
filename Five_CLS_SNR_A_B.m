%%测试M对估计性能的影响
%%整理之后的CLS针对自相关序列
clc;
close all;
clear all;

%% 信号参数
N = 200;
n = 1:1:N;
L = 100;

w0 = 0.4 * pi;   %从pi/25到pi-pi/25每隔pi/25取一个
B =1;
AA = 0:0.02:2; 
phiA = 0.1*pi;
phiB = 0.2*pi;

%% 独立运行
SNR1 = 0;
SNR2 =10;
SNR3 =20;
SNR4 = 40;
T = 1000;
index = 1;
for A = AA 
    wAcu = w0 * ones(1,T);             %真实的信号
    for t = 1:T
        n = 1:1:N;
        sn = (A+B)*exp(- 1j*phiA)*exp(1j*w0*n) +  B*exp(- 1j*phiB)*exp(-1j*w0*n);
        
        xn1 = awgn(sn, SNR1);         %RPHD的含有加性高斯白噪声的信号
        xn2 = awgn(sn, SNR2); 
        xn3 = awgn(sn, SNR3); 
        xn4 = awgn(sn, SNR4); 
        
        for k=1:N-1
            rn1(k)=r(xn1,N,k);
        end
        for k=1:N-1
            rn2(k)=r(xn2,N,k);
        end
        for k=1:N-1
            rn3(k)=r(xn3,N,k);
        end
        for k=1:N-1
            rn4(k)=r(xn4,N,k);
        end
       
        %% CLS1方法的参数
          Rn1 = rn1 - conj(rn1);  %%新构造的序列
          LRn1 = length(Rn1);     %%新构造的序列的长度
         v1_current = zeros(L,1); %x(n+1)
         v1_prev_1 = zeros(L,1);  %x(n)
         v1_prev_2 = zeros(L,1);  %x(n-1)
         Fre_CLS1 = zeros(LRn1,1); %%存储值
         w1 = zeros(LRn1,1);   %%估计的值
         for kk = LRn1:-1:(L+2)

             if kk > (L + 1)
                v1_current = [Rn1(kk:-1:(kk-L+1))];   %x(n+1)
                v1_prev_1 = [Rn1((kk-1):-1:(kk-L))];  %x(n)
                v1_prev_2 = [Rn1((kk-2):-1:(kk-L-1))]; %x(n-1)
                v1_current = v1_current';
                v1_prev_1 = v1_prev_1';
                v1_prev_2 = v1_prev_2';
                w1(kk) = (v1_prev_1'*(v1_current+v1_prev_2))/( (v1_prev_1)'*(v1_prev_1)); 
             end
             Fre_CLS1(kk) = acos(0.5 * w1(kk));
         end
         CLS_gu1(t) = mean(Fre_CLS1(L+2:LRn1));

         %% CLS2的关键参数
         Rn2 = rn2 - conj(rn2);  %%新构造的序列
         LRn2 = length(Rn2);     %%新构造的序列的长度
         v2_current = zeros(L,1); %x(n+1)
         v2_prev_1 = zeros(L,1);  %x(n)
         v2_prev_2 = zeros(L,1);  %x(n-1)
         Fre_CLS2 = zeros(LRn2,1); %%存储值
         w2 = zeros(LRn2,1);   %%估计的值
         for kk = LRn2:-1:(L+2)

             if kk > (L + 1)
                v2_current = [Rn2(kk:-1:(kk-L+1))];   %x(n+1)
                v2_prev_1 = [Rn2((kk-1):-1:(kk-L))];  %x(n)
                v2_prev_2 = [Rn2((kk-2):-1:(kk-L-1))]; %x(n-1)
                v2_current = v2_current';
                v2_prev_1 = v2_prev_1';
                v2_prev_2 = v2_prev_2';
                w2(kk) = (v2_prev_1'*(v2_current+v2_prev_2))/( (v2_prev_1)'*(v2_prev_1)); 
             end
             Fre_CLS2(kk) = acos(0.5 * w2(kk));
         end
         CLS_gu2(t) = mean(Fre_CLS2(L+2:LRn2));
         
         %% CLS3的关键参数
        Rn3 = rn3 - conj(rn3);  %%新构造的序列
        LRn3 = length(Rn3);     %%新构造的序列的长度
        v3_current = zeros(L,1); %x(n+1)
         v3_prev_1 = zeros(L,1);  %x(n)
         v3_prev_2 = zeros(L,1);  %x(n-1)
         Fre_CLS3 = zeros(LRn3,1); %%存储值
         w3 = zeros(LRn3,1);   %%估计的值
         for kk = LRn3:-1:(L+2)

             if kk > (L + 1)
                v3_current = [Rn3(kk:-1:(kk-L+1))];   %x(n+1)
                v3_prev_1 = [Rn3((kk-1):-1:(kk-L))];  %x(n)
                v3_prev_2 = [Rn3((kk-2):-1:(kk-L-1))]; %x(n-1)
                v3_current = v3_current';
                v3_prev_1 = v3_prev_1';
                v3_prev_2 = v3_prev_2';
                w3(kk) = (v3_prev_1'*(v3_current+v3_prev_2))/( (v3_prev_1)'*(v3_prev_1)); 
             end
             Fre_CLS3(kk) = acos(0.5 * w3(kk));
         end
         CLS_gu3(t) = mean(Fre_CLS3(L+2:LRn3));
         
         %% CLS4的关键参数
        Rn4 = rn4 - conj(rn4);  %%新构造的序列
        LRn4 = length(Rn4);     %%新构造的序列的长度
        v4_current = zeros(L,1); %x(n+1)
         v4_prev_1 = zeros(L,1);  %x(n)
         v4_prev_2 = zeros(L,1);  %x(n-1)
         Fre_CLS4 = zeros(LRn4,1); %%存储值
         w4 = zeros(LRn4,1);   %%估计的值
         for kk = LRn4:-1:(L+2)

             if kk > (L + 1)
                v4_current = [Rn4(kk:-1:(kk-L+1))];   %x(n+1)
                v4_prev_1 = [Rn4((kk-1):-1:(kk-L))];  %x(n)
                v4_prev_2 = [Rn4((kk-2):-1:(kk-L-1))]; %x(n-1)
                v4_current = v4_current';
                v4_prev_1 = v4_prev_1';
                v4_prev_2 = v4_prev_2';
                w4(kk) = (v4_prev_1'*(v4_current+v4_prev_2))/( (v4_prev_1)'*(v4_prev_1)); 
             end
             Fre_CLS4(kk) = acos(0.5 * w4(kk));
         end
         CLS_gu4(t) = mean(Fre_CLS4(L+2:LRn4));  
                
    end
    MSE_CLS1(index) = 10*log10(mse(wAcu, CLS_gu1));
    MSE_CLS2(index) = 10*log10(mse(wAcu, CLS_gu2));
    MSE_CLS3(index) = 10*log10(mse(wAcu, CLS_gu3));
    MSE_CLS4(index) = 10*log10(mse(wAcu, CLS_gu4));

    index
    index = index + 1;
end
plot(AA, MSE_CLS1, '-', 'LineWidth', 2);
hold on;
plot(AA, MSE_CLS2, '-', 'LineWidth', 2);
plot(AA, MSE_CLS3, '-', 'LineWidth', 2);
plot(AA, MSE_CLS4, '-', 'LineWidth', 2);
%仿真格式部分
legend('SNR=0dB','SNR=10dB','SNR=20dB','SNR=40dB');%用指定的文字在当前坐标轴中对所给数据的每一部分显示一个图例
xlabel('\fontname{Times New Roman}A-B', 'FontWeight','bold');%字体Times New Roman，加粗
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;