function [x_est] = sft(x,k)
% 输入序列以及稀疏度，返回sft结果
n = length(x);
B = 2^round(log(round(sqrt(n*k)))/log(2));%使得B是2的幂次，保证整除
delta = 1/n;%暂定
d = 2;
L = 5*round(log(n)/log(2));
theta = 2/3; % 自己添加的参数，超过theta*L次后认为存在这一值
w = min(n,round(B*log(n/delta)));

[g_window,~,g_window_f] = gauss_window(n,1/(2*B),1/B,delta,min(n,round(B*log(n/delta))));

loc_loop = zeros(1,n);

% 先进行定位循环
for count = 1:1:L
    % step1:permute,将大频点打乱
    sigma = randi([0,n/2-1])*2+1;%必须是奇数以保证有逆
    tau = randi([0,n-1]);
    x_p = zeros(size(x));
    for ii = [1:floor(w/2),n-floor(w/2):n]
        x_p(ii) = x(mod(sigma*(ii-1)+tau,n)+1);
    end
    
    % step2:卷上频域窗函数
    x_g = zeros(size(x));
    x_g([1:floor(w/2),n-floor(w/2):end]) = g_window([1:floor(w/2),n-floor(w/2):end]).*x_p([1:floor(w/2),n-floor(w/2):end]);
    
    % step3:频域降采样
    y = zeros(1,B);
    for ii = 1:1:B
        y(ii) = sum(x_g(ii:B:end));
    end
    z = fft(y);
    
    % 找到幅度最大的d*k个
    [~,J] = maxk(abs(z),d*k);
    for ii = 1:1:length(J)
        list = inv_hash(n,sigma,B,J(ii)-1)+1;
        list = list(list>0);
        loc_loop(list) = loc_loop(list)+1;
    end
end
loc_loop(loc_loop<theta*L) = 0;
loc_loop(loc_loop>=theta*L) = 1;

list = find(loc_loop);
est_loop = zeros(length(list),L);

% 再进行值判决，多次循环取中间值！（去除频点过于接近导致的错误值）
for count = 1:1:L
    % step1:permute,将大频点打乱
    sigma = randi([0,n/2-1])*2+1;%必须是奇数以保证有逆
    tau = randi([0,n-1]);
    x_p = zeros(size(x));
    for ii = [1:floor(w/2),n-floor(w/2):n]
        x_p(ii) = x(mod(sigma*(ii-1)+tau,n)+1);
    end
    
    % step2:卷上频域窗函数
    x_g = zeros(size(x));
    x_g([1:floor(w/2),n-floor(w/2):end]) = n*g_window([1:floor(w/2),n-floor(w/2):end]).*x_p([1:floor(w/2),n-floor(w/2):end]);
    
    % step3:频域降采样
    y = zeros(1,B);
    for ii = 1:1:B
        y(ii) = sum(x_g(ii:B:end));
    end
    z = fft(y);
    
    % 估计大小
    est_loop(:,count) = z(hash(n,sigma,B,list-1)+1).*exp(-1i*2*pi*tau.*(list-1)/n)./g_window_f(1+mod(mod(sigma*(list-1),n)-round((hash(n,sigma,B,list-1))*(n/B)),n));
end
est_loop= sort(est_loop,2);
x_est = zeros(1,n);
x_est(list) = est_loop(:,round(L/2));
end