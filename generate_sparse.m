function [x,t,k_list] = generate_sparse(n,k,amin,amax,anoise)
% amin,amax 表示频域幅度随机值的范围,anoise为其它点的噪声幅度
if(~exist('amax','var'))&&(~exist('amin','var'))
    amax = 1;  % 如果未出现该变量，则对其进行赋值
    amin = 0.5;
end
if(~exist('anoise','var'))
    anoise = 0;
end
if k>n
    x = [];
else
    k_list = randsample(n,k);
    t = 0:1:n-1;
    x_f = rand(size(t))*anoise;
    for ii = 1:1:length(k_list)
        x_f(k_list(ii)) = (rand*(amax-amin)+amin)*exp(1i*rand*2*pi);
    end
    x = ifft(x_f);
end