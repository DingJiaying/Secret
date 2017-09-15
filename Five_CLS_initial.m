%%����֮���CLS
clc;
close all;
clear all;

%% �źŲ���
N = 200;  %���г���
n = 1:1:N; 
L = 100;      %������
w0 = 0.02*pi;


A =1;
B = 2;
phiA = 0.1*pi;
phiB = 0.2*pi;
sn = A*exp(- 1j*phiA)*exp(1j*w0*n) +  B*exp(- 1j*phiB)*exp(-1j*w0*n);


%% ��������
SNR = 20;
T = 100;     %%�������еĴ���
wAcu = w0 * ones(1,T);             %��ʵ���ź�
for t = 1:T
    xn = awgn(sn,SNR);
    for k=1:N-1
        rn(k)=r(xn,N,k);
    end
    Rn = rn - conj(rn);  %%�¹��������
    LRn = length(Rn);  %%�¹�������еĳ���
    %% CLS�����Ĳ���
     v_current = zeros(L,1); %x(n+1)
     v_prev_1 = zeros(L,1);  %x(n)
     v_prev_2 = zeros(L,1);  %x(n-1)
     Fre_CLS = zeros(LRn,1); %%�洢ֵ
     w = zeros(LRn,1);   %%���Ƶ�ֵ
     for kk = LRn:-1:(L+2)
      
         if kk > (L + 1)
            v_current = [Rn(kk:-1:(kk-L+1))];   %x(n+1)
            v_prev_1 = [Rn((kk-1):-1:(kk-L))];  %x(n)
            v_prev_2 = [Rn((kk-2):-1:(kk-L-1))]; %x(n-1)
            v_current = v_current';
            v_prev_1 = v_prev_1';
            v_prev_2 = v_prev_2';
%             w(kk) = (v_prev_1'*(v_current+v_prev_2))/( norm(v_prev_1).^2); 
            w(kk) = (v_prev_1'*(v_current+v_prev_2))/( (v_prev_1)'*(v_prev_1)); 
         end
         Fre_CLS(kk) = acos(0.5 * w(kk));
     end
     CLS_gu(t) = mean(Fre_CLS(L+2:LRn));
      
end
MSE_CLS = 10*log10(mse(wAcu, CLS_gu))
 