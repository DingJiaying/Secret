%%����M�Թ������ܵ�Ӱ��
%%����֮���CRPHD������������
clc;
close all;
clear all;

%% �źŲ���
L = 100;

w0 = 0.3*pi;
A =1;
B = 2;
phiA = 0.1*pi;
phiB = 0.2*pi;

SNR1 = 0;
SNR2 =10;
SNR3 =20;
SNR4 = 40;

%% ��������
T = 1000;
NN = 100:20:600;

index = 1;
for N = NN
     wAcu = w0 * ones(1,T);             %��ʵ���ź�
    
    n = 1:1:N;
    sn = A*exp(- 1j*phiA)*exp(1j*w0*n) +  B*exp(- 1j*phiB)*exp(-1j*w0*n);
    for t = 1:T
        xn1 = awgn(sn, SNR1);         %RPHD�ĺ��м��Ը�˹���������ź�
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
        
        %% SCRPHD1�����Ĳ��� 
         Rn1 = rn1 - conj(rn1);  %%�¹��������
         LRn1 = length(Rn1);     %%�¹�������еĳ���
         Fre_SCRPHD1 = zeros(LRn1,1); %%�洢ֵ
         for kk = LRn1:-1:(L+2)      
             if kk > (L + 1)            
                SCRPHD_A1 = real(AR(Rn1(kk:-1:(kk-L+1)), L));
                SCRPHD_B1 = real(BR(Rn1(kk:-1:(kk-L+1)), L));
                argument_r1 = (SCRPHD_B1+sqrt(SCRPHD_B1^2+8*SCRPHD_A1^2))/(4*SCRPHD_A1);  
                if (argument_r1>1)
                    argument_r1 = 1;
                end
                if (argument_r1<-1)
                   argument_r1 = -1;
                end
             end
             Fre_SCRPHD1(kk) = acos(argument_r1);
         end
         SCRPHD_gu1(t) = mean(Fre_SCRPHD1(L+2:LRn1));
    
    %% SCRPHD2�����Ĳ���  
        Rn2 = rn2 - conj(rn2);  %%�¹��������
        LRn2 = length(Rn2);     %%�¹�������еĳ���
         Fre_SCRPHD2 = zeros(LRn2,1); %%�洢ֵ
         for kk = LRn2:-1:(L+2)      
             if kk > (L + 1)            
                SCRPHD_A2 = real(AR(Rn2(kk:-1:(kk-L+1)), L));
                SCRPHD_B2 = real(BR(Rn2(kk:-1:(kk-L+1)), L));
                argument_r2 = (SCRPHD_B2+sqrt(SCRPHD_B2^2+8*SCRPHD_A2^2))/(4*SCRPHD_A2);  
                if (argument_r2>1)
                    argument_r2 = 1;
                end
                if (argument_r2<-1)
                   argument_r2 = -1;
                end
             end
             Fre_SCRPHD2(kk) = acos(argument_r2);
         end
         SCRPHD_gu2(t) = mean(Fre_SCRPHD2(L+2:LRn2));
        
        %% SCRPHD3�����Ĳ��� 
         Rn3 = rn3 - conj(rn3);  %%�¹��������
         LRn3 = length(Rn3);     %%�¹�������еĳ���
         Fre_SCRPHD3 = zeros(LRn3,1); %%�洢ֵ
         for kk = LRn3:-1:(L+2)      
             if kk > (L + 1)            
                SCRPHD_A3 = real(AR(Rn3(kk:-1:(kk-L+1)), L));
                SCRPHD_B3 = real(BR(Rn3(kk:-1:(kk-L+1)), L));
                argument_r3 = (SCRPHD_B3+sqrt(SCRPHD_B3^2+8*SCRPHD_A3^2))/(4*SCRPHD_A3);  
                if (argument_r3>1)
                    argument_r3 = 1;
                end
                if (argument_r3<-1)
                   argument_r3 = -1;
                end
             end
             Fre_SCRPHD3(kk) = acos(argument_r3);
         end
         SCRPHD_gu3(t) = mean(Fre_SCRPHD3(L+2:LRn3));
         
        %% SCRPHD4�����Ĳ���
        Rn4 = rn4 - conj(rn4);  %%�¹��������
        LRn4 = length(Rn4);     %%�¹�������еĳ���
         Fre_SCRPHD4 = zeros(LRn4,1); %%�洢ֵ
         for kk = LRn4:-1:(L+2)      
             if kk > (L + 1)            
                SCRPHD_A4 = real(AR(Rn4(kk:-1:(kk-L+1)), L));
                SCRPHD_B4 = real(BR(Rn4(kk:-1:(kk-L+1)), L));
                argument_r4 = (SCRPHD_B4+sqrt(SCRPHD_B4^2+8*SCRPHD_A4^2))/(4*SCRPHD_A4);  
                if (argument_r4>1)
                    argument_r4 = 1;
                end
                if (argument_r4<-1)
                   argument_r4 = -1;
                end
             end
             Fre_SCRPHD4(kk) = acos(argument_r4);
         end
         SCRPHD_gu4(t) = mean(Fre_SCRPHD4(L+2:LRn4));
    end
    MSE_SCRPHD1(index) = 10*log10(mse(wAcu, SCRPHD_gu1));
    MSE_SCRPHD2(index) = 10*log10(mse(wAcu, SCRPHD_gu2));
    MSE_SCRPHD3(index) = 10*log10(mse(wAcu, SCRPHD_gu3));
    MSE_SCRPHD4(index) = 10*log10(mse(wAcu, SCRPHD_gu4));    
    index
    index = index + 1;
end
plot(NN, MSE_SCRPHD1, '-^', 'LineWidth', 2);
hold on;
plot(NN, MSE_SCRPHD2, '-o', 'LineWidth', 2);
plot(NN, MSE_SCRPHD3, '-*', 'LineWidth', 2);
plot(NN, MSE_SCRPHD4, '-d', 'LineWidth', 2);
%�����ʽ����
legend('SNR=0dB','SNR=10dB','SNR=20dB','SNR=40dB');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}N', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;