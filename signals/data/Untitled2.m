

w0 = 1; %hz 
w1 = 100; %hzs
r = 2*pi / 1000;
x = -pi:r:pi; 
Y =  0.5*(1 + cos(w0*x).*cos(w1*x));
Y = Y.';
%Y(1)
%A = [x(:), Y(:) ];
%A
save(sprintf("C:\\Users\\Owner\\Desktop\\michael_david\\yvec%d_%d.txt", w0, w1), 'Y', '-ascii')
