function [dI0dx,dI0dy]=gradofI(I)
b0=gpuArray(single([1/120 13/60 11/20 13/60 1/120]));
b1=gpuArray(single([-1/24 -5/12 0 5/12 1/24]));
dI0dx=conv2(fdeconv(padarray(I,[3 3],'replicate'),b0),b1);
dI0dy=conv2(fdeconv(padarray(I,[3 3],'replicate')',b0),b1)';

dI0dx=dI0dx(4:end-3,4:end-3);
dI0dy=dI0dy(4:end-3,4:end-3);

dI0dx=-dI0dx;
dI0dy=-dI0dy;
