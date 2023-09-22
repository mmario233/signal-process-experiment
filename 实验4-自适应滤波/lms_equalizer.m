function [ z ] = lms_equalizer( y,info,delta )
%LMS�㷨����Ӧ�˲���ʵ��
estimated_c=[0 0 0 0 0 1 0 0 0 0 0]; %��ʼ��ͷϵ��(����Ӧ���Ǻ�channel���ŵ�������ͬ=11)
K=5;                                 %K=��length��estimated_c��-1��/2
for k=1:size(y,2)-2*K              %channel�з��ز���len�ĳ���Ҳ��5������K��ѡ����ǻ���len����ҪK=len
     y_k=y(k:k+2*K);                 %��ȡ��Ԫ��һ��11��
     z_k=estimated_c*y_k';           %����ͷϵ������Ԫ��˺����
     e_k=info(k)-z_k;                %������
     estimated_c=estimated_c+delta * e_k * y_k;%����У����ͷϵ��
     z(k)=z_k;                       %������������Ԫ����
end
%���e=d-y����yָ�����ŵ��������źţ����������ź�dʹ�õ��������ź�info
%�Ƚ������ʶ���ֻ�Ƚ�info�ĳ���N
%size(y,2)����y������=N+L-1��L=2len+1��,��size(y,2)-2*K=N+L-(2K+1)=N+2(len-K)
%���ｫ���и��׶ˣ����len>K����k=1��N+2(len-K)�ᳬ��info�ĳ���N.
end