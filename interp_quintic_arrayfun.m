function J=interp_quintic_arrayfun(F,X,Y)
    h=size(F,1);
    w=size(F,2);
    x=max(min(X,w),1);
    y=max(min(Y,h),1);
    x_int=floor(x);
    y_int=floor(y);
    dx=(x-x_int);
    dy=(y-y_int);

    J=arrayfun(@interp_point,x_int,y_int,dx,dy);

    function f=interp_point(x,y,dx,dy)
        f=F(y,x,6,1)+dx*(F(y,x,6,2)+dx*(F(y,x,6,3)+dx*(F(y,x,6,4)+dx*(F(y,x,6,5)+dx*F(y,x,6,6))))); 
        f=f*dy+F(y,x,5,1)+dx*(F(y,x,5,2)+dx*(F(y,x,5,3)+dx*(F(y,x,5,4)+dx*(F(y,x,5,5)+dx*F(y,x,5,6)))));
        f=f*dy+F(y,x,4,1)+dx*(F(y,x,4,2)+dx*(F(y,x,4,3)+dx*(F(y,x,4,4)+dx*(F(y,x,4,5)+dx*F(y,x,4,6)))));
        f=f*dy+F(y,x,3,1)+dx*(F(y,x,3,2)+dx*(F(y,x,3,3)+dx*(F(y,x,3,4)+dx*(F(y,x,3,5)+dx*F(y,x,3,6)))));
        f=f*dy+F(y,x,2,1)+dx*(F(y,x,2,2)+dx*(F(y,x,2,3)+dx*(F(y,x,2,4)+dx*(F(y,x,2,5)+dx*F(y,x,2,6)))));
        f=f*dy+F(y,x,1,1)+dx*(F(y,x,1,2)+dx*(F(y,x,1,3)+dx*(F(y,x,1,4)+dx*(F(y,x,1,5)+dx*F(y,x,1,6)))));
    end
end