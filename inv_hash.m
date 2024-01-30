function i = inv_hash(n,sigma,B,h)
% 寻找i使得sigma*i = h*(n/B) mod n
% 一共有最多floor(n/B)+1个点
% logn复杂度

% 首先寻找sigma的modn逆元
% sigma*inv + x*n = 1
% 扩展欧几里得算法
a = sigma;
b = n;
inv = 1;%a中有几个sigma
x = 0;% a中有几个n
n_n = 1; % b中有几个n
n_sig = 0;% b中有几个sigma
while a>1
    % 先找a中的b
    if b>1
        inv = inv-n_sig*floor(a/b);
        x = x-n_n*floor(a/b);
        a = mod(a,b);
    else
        inv = inv-(a-1)*n_sig;
        break
    end

    % 再反过来
    n_n = n_n-x*floor(b/a);
    n_sig = n_sig-inv*floor(b/a);
    b = mod(b,a);
end
inv = mod(inv,n);
i = ones(1,floor(n/B)+1)*(-1);
si = ceil((h-1/2)*n/B);
count = 1;
while si<ceil((h+1/2)*n/B)
    i(count) = mod(inv*si,n);
    si = si+1;
    count = count+1;
end
end