function [ info ] = random_binary( N )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if nargin == 0      %nargin��ʾ�����õĺ�������������ĸ���
    N=10000;         %���û�������������ָ����Ϣ����Ϊ10000����Ԫ
end
for i=1:N
    temp=rand;
    if(temp<0.5)
        info(i)=-1;  %1/2�ĸ���
    else
        info(i)=1;
    end
end
end