%%����M�Թ������ܵ�Ӱ��
%%����֮���CRPHD������������
clc;
close all;
clear all;

%% �źŲ���
N = 200;
n = 1:1:N;
L = 100;

A =1;
B = 2;
phiA = 0.1*pi;
phiB = 0.2*pi;
w0 = 0.2*pi;   %��pi/25��pi-pi/25ÿ��pi/25ȡһ��

sn = A*exp(- 1j*phiA)*exp(1j*w0*n) +  B*exp(- 1j*phiB)*exp(-1j*w0*n);

LL = 10:10:N-1;

%% ��������
T = 1000;
SNR = 20;
index = 1;
for L = LL
    wAcu = w0 * ones(1,T);             %��ʵ���ź�
    for t = 1:T
        xn = awgn(sn,SNR);
        for k=1:N-1
            rn(k)=r(xn,N,k);
        end
        Rn = rn - conj(rn);  %%�¹��������
        LRn = length(Rn);  %%�¹�������еĳ���
    %% SCRPHD�����Ĳ���    
     Fre_SCRPHD = zeros(LRn,1); %%�洢ֵ
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

%�����ʽ����
legend('SCRPHD');%��ָ���������ڵ�ǰ�������ж��������ݵ�ÿһ������ʾһ��ͼ��
xlabel('\fontname{Times New Roman}L', 'FontWeight','bold');%����Times New Roman���Ӵ�
ylabel('\fontname{Times New Roman}MSE (dB)', 'FontWeight','bold');
grid on;
hold off;