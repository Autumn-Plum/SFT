function h = hash(n,sigma,B,i)
% 如果输入的是序号，需要先减一！
h = round(mod(i*sigma,n)*(B/n));
h(h==B) = 0;
end