k=[190.33 43.77 9.73 60.54 49.01 9.04;135.2 36.4 10.47 44.16 36.49 3.94;95.21 22.83 9.3 22.44 22.81 2.8;104.78 25.11 6.4 9.89 18.17 3.25;128.41 27.63 8.94 12.58 23.99 3.27;145.68 32.83 17.79 27.29 39.09 3.47;159.37 33.38 18.37 11.81 25.29 5.22; 116.22 29.57 13.24 13.76 21.75 6.04;221.11 38.64 12.53 115.65 50.82 5.89;144.98 29.12 11.67 42.6 27.3 5.74;169.92 32.75 12.72 47.12 34.35 5;153.11 23.09 15.62 23.54 18.18 6.39;144.92 21.26 16.96 19.52 21.75 6.73;140.55 21.5 17.64 19.19 15.97 4.94;115.84 30.26 12.2 33.61 33.77 3.85;101.18 23.26 8.46 20.2 20.5 4.3];
F=k/sum(sum(k));
Dn=diag(sum(F,2));
Dp=diag(sum(F,1));
F=F-Dn*ones(size(k,1),1)*(Dp*ones(size(k,2),1))';
S=F'*Dn^-1*F*Dp^-1;
Q=F*Dp^-1*F'*Dn^-1;
[VS,DS]=eig(S);
[VQ,DQ]=eig(Q);
%normalize the principal axis, calculate the factors
u=VS*(VS'*Dp^-1*VS)^-0.5;
r=Dp^-1*u;
v=VQ*(VQ'*Dn^-1*VQ)^-0.5;
for i=1:min(size(k,1),size(k,2))
    if F'*Dn^-1*v(:,i)./u(:,i)<0
        v(:,i)=v(:,i)*-1;
    end
end
s=Dn^-1*v;
for i=1:min(size(k,1),size(k,2))
    if F'*Dn^-1*v(:,i)./u(:,i)<0
        error('the direction of u and v is not corresponding!');
    elseif Dp^-1*F'*s(:,i)./r(:,i)<0
        error('the direction of r and s is not corresponding!');
    end
end
r2=r(:,1:2)*DS(1:2,1:2)^0.5;
s2=s(:,1:2)*DS(1:2,1:2)^0.5;
%plotting
plot(r2(:,1),r2(:,2),'r.',s2(:,1),s2(:,2),'b.')
num1={'1','2','3','4','5','6'};
num2={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'};
text(r2(:,1)+0.01,r2(:,2)+0.01,num1);
text(s2(:,1)+0.01,s2(:,2)+0.01,num2);
%X4, the cost of housing, is the first kind of expenditure, the people live in the provinces with large population spend more on housing than those living in provinces with small population.
%X1, X3, X6 are the second kind of expenditrue, people live in the provinces which is more developed or in the south-eastern part of China spend more on food, fuel and culture than other people.
%The last type of expenditure includes X2 and X5, clothing and daily use.
%kk=k([9,1,2,10,11,15,3,16,6,12,13,5,14,4,8,7],[4,5,2,1,6,3]);
