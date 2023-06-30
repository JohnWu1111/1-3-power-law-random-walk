clear;
% close all;
clc;
format long
tic;

t = 1:100;
nt = length(t);
value = zeros(nt,1);
for i = 1:nt
    value(i) = SimpInt(0, 1, 1e4, t(i));
end

figure
plot(t, value)
% plot(log(t), log(value))
plot(t, log(value))

b = 1e-3:1e-3:1;
a = -log(b)./b;
figure
plot(log(a),log(b))

fitx = log(t);
fity = log(value);

toc;

function S = SimpInt(a, b, num, t)
n = 2*num+1;
h = (b-a)/(n-1);
x = a+h:h:b;
n = n-1;
y = f(x, t);
S = (h/3)*(y(1) + 4*sum(y(2:2:n-1)) + 2*sum(y(3:2:n-2)) + y(n));
end

function y = f(x, t)
y = exp(-t*((-log(x)).^(0)))./x;
end