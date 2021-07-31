
function [x, Y, mu, sig, w0, mu2, sig2,  w1] = create_signal(mu, sig, w0, mu2, sig2,  w1)
    r = 2*pi / 10000;
    x = -pi:r:pi; 
    Y =   gaussmf(x,[sig2, mu2]) +  gaussmf(x,[sig, mu]).*exp(1i*w0*x);
    Y = 0.5 *( 1+ Y.*exp(1i*w1*x));
    Y = Y.';
end