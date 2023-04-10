clear all
close all
clc

% m = 1;
% n = 2;
% p = 2;
% q = 2;

N = 1e4;
U = 30;
L = 0;
x = linspace(L, U, N);

ln = log(1 + x);

A = [1,1];
B = [];
C = [1];
D = [0];

% for i = 1:N
%     meijer(i) = meijerG(A, B, C, D, x(i))
% end

an = [1, 1];
An = [1, 1];
ap = [];
Ap = [];
bm = [1];
Bm = [1];
bq = [0];
Bq = [1];

fox = zeros(1,N);
for i = 1:N
    fox(i) = HFox(an, An, ap, Ap, bm, Bm ,bq, Bq, x(i));
end

% plot(x, ln, 'r--', x, fox, 'b:')

err = abs(ln - fox);
plot(x,err)
