function [x,t,x_f] = gauss_window(n,eps1,eps2,delta,w)
%参数按照文章里的定义,n为2的幂次
%认为卷积的矩形长度为(eps1+eps2)*n(宽于平坦宽度！要注意),高斯的半宽度为(eps1-eps2)*n/2
f = 0:1:n-1;
f(n/2+1:n) = f(n/2+1:n)-n;
gauss = exp(-log(1/delta)*f.^2./((eps1-eps2)*n/2)^2);%频域高斯
window = zeros(size(f));
window(-round((eps1+eps2)*n/2)<f&f<round((eps1+eps2)*n/2)) = 1;
x_f = ifft(fft(window).*fft(gauss));
x_f = x_f(1:n);
x_f = x_f/max(x_f);
x_f(abs(x_f)<1e-15) = 0;
x = real(ifft(x_f));
x(-round(w/2)>f|f>round(w/2)) = 0;
t = 0:1:n-1;
end