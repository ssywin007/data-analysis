k=[58 19 1 33 110 163;5 1 1 1 31 25;20 3 2 10 78 121;1 0 0 0 12 14;1 0 0 2 28 3;1 0 1 4 16 10;6 2 1 7 19 17;3 0 0 10 19 12;3 0 1 0 30 31];
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
num2={'1','2','3','4','5','6','7','8','9'};
text(r2(:,1)+0.01,r2(:,2)+0.01,num1);
text(s2(:,1)+0.01,s2(:,2)+0.01,num2);
%the outcoming is the same as P64
%kk=k([1,7,3,8,4,2,9,6,5],[2,1,4,6,5,3]);
