clear,clc
directory='images/';
files=dir([directory '*.tiff']);

I0=imread([directory files(1).name]);
I0=gpuArray(single(rgb2gray(I0)));

% Locations of centers
xc=[628 1264 226 1395]';
yc=[139 191 666 886]';

hws=16; % Half window size

c=sub2ind(size(I0)-2*hws,yc-hws,xc-hws);

image(I0)
colormap(gray(256))
hold on
plot(xc(:)'+[-hws hws hws -hws -hws]',yc(:)'+[-hws -hws hws hws -hws]')
hold off
drawnow
%%
queuesize=20000;
I0_info=calculate_I0_parameters(I0,hws,xc,yc,queuesize);
fprintf('Image %s loaded.\n',files(1).name)
%%
queuesize=20000;
Mp=repmat(permute(single([1 0 0;0 1 0]),[3 4 1 2]),[size(xc) 1 1]);
u=zeros(numel(xc),numel(files));
v=u;
CC=u+1;
for i=2:numel(files)
    I1=imread([directory files(i).name]);
    I1=gpuArray(single(rgb2gray(I1)));
%     [u,v,M]=int_disp(I0,I1,hws,10,10,100000);
%     Mp(:,:,1,3)=u(c);
%     Mp(:,:,2,3)=v(c);
    [Mp,cc]=IC_method(I0_info,I1,Mp,hws,xc,yc,10,1e-4,queuesize);
    CC(:,i)=gather(cc);
    u(:,i)=Mp(:,1,1,3);
    v(:,i)=Mp(:,1,2,3);
    if mod(i,100)==0
        plot(u(:,1:i)')
        drawnow
    end
end
hist(counterout(:),0:20)