function Mp2=inverseComp2(Mp,d)
Mp2=Mp;
Mp2(:,1,1)=1+d(:,3);
Mp2(:,1,2)=d(:,4);
Mp2(:,1,3)=d(:,1);
Mp2(:,2,1)=d(:,5);
Mp2(:,2,2)=1+d(:,6);
Mp2(:,2,3)=d(:,2);
Mp2(:,3,3)=1;
Mp(:,3,3)=1;

Mp2=permute(pagefun(@mrdivide,permute(Mp,[2 3 1]),permute(Mp2,[2 3 1])),[3 1 2]);
Mp2=Mp2(:,1:2,:);