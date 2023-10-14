function F=quintic_coeffs(I)
h=gpuArray(single([1/120 13/60 11/20 13/60 1/120]));
F=zeros([size(I) 6 6],'single','gpuArray');
Ipad=padarray(I,[5,5],'pre','replicate');
Ipad=padarray(Ipad,[6,6],'post','replicate');
Idec=fdeconv(Ipad,h)';
Idec=fdeconv(Idec,h)';
QK=[1/120 13/60 11/20 13/60 1/120 0;
    -1/24 -5/12 0 5/12 1/24 0;
    1/12 1/6 -1/2 1/6 1/12 0;
    -1/12 1/6 0 -1/6 1/12 0;
    1/24 -1/6 1/4 -1/6 1/24 0;
    -1/120 1/24 -1/12 1/12 -1/24 1/120];
for i=1:6
    for j=1:6
        Iconv=conv2(fliplr(QK(i,:)),fliplr(QK(j,:)),(Idec),'valid');
        F(:,:,i,j)=Iconv(2:end-1,2:end-1);
    end
end