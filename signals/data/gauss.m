

w0 = 1;
w1 = 20;

[x0, Y0, mu, sig, w0 , mu2, sig2, w2]  = create_signal(0, 4, w0, 0, 4, w1);
for i = 1:20:100
    [x, Y, mu, sig, w0 , mu2, sig2, w2]  = create_signal(0, 4, i * w0, 0, 4, i*w1);
    Y0 = Y0 + Y;
end 
save(sprintf("C:\\Users\\Owner\\Desktop\\michael_david\\Gauss.txt"), 'Y0', '-ascii');
fs = 10000;
plot(x0,Y0);

X = fft(Y0);
X = fftshift(abs(X));
n = length(X);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
%powershift = abs(X).^2/n;     % zero-centered power
%plot(fshift,X);