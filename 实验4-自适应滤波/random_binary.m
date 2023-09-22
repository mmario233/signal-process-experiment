function [ info ] = random_binary( N )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin == 0      %nargin表示所引用的函数的输入参数的个数
    N=10000;         %如果没有输入参数，则指定信息序列为10000个码元
end
for i=1:N
    temp=rand;
    if(temp<0.5)
        info(i)=-1;  %1/2的概率
    else
        info(i)=1;
    end
end
end