clear;
clc;
rng(2024);
%运行时间较长
n_list = 2.^[10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26];
t_run = zeros(size(n_list));

for n = n_list
%n = 2^15;
k = 25;
B = 2^round(log(floor(sqrt(n*k)))/log(2));%使得B是2的幂次，保证整除
delta = 1/n;%暂定
d = 2;
L = round(log(n)/log(2));
theta = 2/3; % 自己添加的参数，超过theta*L次后认为存在这一值

w = min(n,round(B*log(n/delta)));

[x,t,k_list] = generate_sparse(n,k,0.5,1,0.01);
x_f = fft(x);

tic;
x_est = sft(x,k);

t_run(log2(n)-9) = toc;
n
t_run(log2(n)-9)
end

figure;
loglog(n_list,t_run,"*-","linewidth",2)
hold on;
loglog(n_list,t_run(end-5)/n_list(end-5)*n_list,"linewidth",2)
legend(["运算时间","O(n)曲线"],"location","northwest");
xlabel("序列长度")

figure;
semilogx(n_list,t_run,"*-","linewidth",2)
hold on;
semilogx(n_list,t_run(end-5)/n_list(end-5)*n_list,"linewidth",2)
legend(["运算时间","O(n)曲线"],"location","northwest");
xlabel("序列长度")

