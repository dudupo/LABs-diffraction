
fs = 1000;


hold on
% cm1hz1khz
A = { Z, Y, U, M, Q};
for S = 1:5
X = fft(table2array(A(S)));
X = fftshift(abs(X));
n = length(X);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(X).^2/n;     % zero-centered power
plot(fshift,powershift);
for c = 1:n
    if powershift(c)> 0.5*(10^(-3))
       disp(fshift(c))
    end
end
end
hold off
shg;