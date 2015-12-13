%%import the pictures, print the original pictures

DW=imread('G:\soon\大四下\数据分析基础\DW.jpg');
DW=DW(:,:,1);
carpet=imread('G:\soon\大四下\数据分析基础\carpet.jpg');
carpet=carpet(:,:,1);
lion=imread('G:\soon\大四下\数据分析基础\lion.jpg');
lion=lion(:,:,1);
imshow(DW);
figure;imshow(carpet);
figure;imshow(lion);

%%mix the pictures random
%generate the coefficiencies
index1=rand(3,1);
index1=index1/norm(index1);
index2=rand(3,1);
index2=index2/norm(index2);
index3=rand(3,1);
index3=index3/norm(index3);
index=[index1';index2';index3';];
%get original signal 's' and mixed signal 'x'
DW1=reshape(DW,size(DW,1)*size(DW,2),1);
carpet1=reshape(carpet,size(carpet,1)*size(carpet,2),1);
lion1=reshape(lion,size(lion,1)*size(lion,2),1);
s=single([DW1';carpet1';lion1']);
x=index*s;
%print mixed pictures
mix1=x(1,:)';
mix1=uint8(mix1);
mix1=reshape(mix1,size(DW,1),size(DW,2));
figure; imshow(mix1);
mix2=x(2,:)';
mix2=uint8(mix2);
mix2=reshape(mix2,size(carpet,1),size(carpet,2));
figure; imshow(mix2);
mix3=x(3,:)';
mix3=uint8(mix3);
mix3=reshape(mix3,size(lion,1),size(lion,2));
figure; imshow(mix3);

%%solve the original signal by maxmizing the L
B=rand(3,3);
for i=1:3
    B(i,:)=B(i,:)/norm(B(i,:));
end
if mean(mean(B*x))>50
    L=log(abs(det(B)))+mean(sum((-2)*(B*x-log(2))));         %%to prevent the lnf
    %L=log(abs(det(B)))+mean(sum(((B*x)-log(2))-(B*x).*(B*x)/2));
else
    L=log(abs(det(B)))+mean(sum((-2)*log(cosh(B*x))));
    %L=log(abs(det(B)))+mean(sum(log(cosh(B*x))-(B*x).*(B*x)/2));
end
L0=L-1;
n=0;
while L-L0>0.01 && n<1000
    gradB=(eye(3)+(-2)*(tanh(B*x)*x'*B')/size(x,2))*B;
    %gradB=(eye(3)+((tanh(B*x)-B*x)*x'*B')/size(x,2))*B;
    deltaB=0.01*(gradB>0)+(-0.01)*(gradB<0);                       %the mutation of B
    B=B+deltaB;
    L0=L;
    if mean(mean(B*x))>50
        L=log(abs(det(B)))+mean(sum((-2)*(B*x-log(2))));
        %L=log(abs(det(B)))+mean(sum(((B*x)-log(2))-(B*x).*(B*x)/2));
    else
        L=log(abs(det(B)))+mean(sum((-2)*log(cosh(B*x))));
        %L=log(abs(det(B)))+mean(sum(log(cosh(B*x))-(B*x).*(B*x)/2));
    end
    n=n+1;
end
%print the solved original pictures
B=B-deltaB;
%normalize the B
for i=1:3
    B(1,:)=B(1,:)/norm(B(1,:));
end
s1=B*x;
if sum(s1(1,:)>0)>sum(s1(1,:)<0)
    s11=reshape(uint8(s1(1,:)),size(DW,1),size(DW,2));
else
    s11=reshape(uint8(-s1(1,:)),size(DW,1),size(DW,2));           %we may get the opposite signals
end
figure; imshow(s11);
if sum(s1(2,:)>0)>sum(s1(2,:)<0)
    s12=reshape(uint8(s1(2,:)),size(DW,1),size(DW,2));
else
    s12=reshape(uint8(-s1(2,:)),size(DW,1),size(DW,2));
end
figure; imshow(s12);
if sum(s1(3,:)>0)>sum(s1(3,:)<0)
    s13=reshape(uint8(s1(3,:)),size(DW,1),size(DW,2));
else
    s13=reshape(uint8(-s1(3,:)),size(DW,1),size(DW,2));
end
figure; imshow(s13);

%%solve the original signal by FastICA Method1
%centralize and whiten the mixed signals
x1=x-mean(x,2)*ones(1,size(x,2));
[V,D]=eig(x1*x1'/size(x,2));
x1=V*D^(-1/2)*V'*x1;

w=zeros(3,3);
pos=1;  
while min(mean(w,2))==0
    for i = pos:3
        wnew=rand(1,3);
        wold=wnew-[1,1,1];
        cnum=0;
        while norm(wold-wnew)>0.01 || norm(wold+wnew)>0.01
            wold=wnew;
            wnew=sum((x1.*(ones(3,1)*tanh(wnew*x1))),2)'/size(x1,2)-sum(sech(wnew*x1).*sech(wnew*x1),2)/size(x1,2)*wnew;
            wnew=wnew/norm(wnew);
            cnum=cnum+1;
            if cnum>1000
                break
            end
        end
        for j=1:i
            if abs(wnew*w(j,:)'-1)<0.01 || abs(wnew*w(j,:)'+1)<0.01          %judge whether wnew has already existed
                break
            end
        end
        if j==i
            wminus=[0,0,0];
            for k=1:i-1 
                wminus=wminus+wnew*w(k,:)'*w(k,:);               %orthogonalize the w
            end
            w(i,:)=wnew-wminus;
            w(i,:)=w(i,:)/norm(w(i,:));
            pos=i+1;                          %record the position of solved w
        end
    end
end
%print the solved original signals
a=V*D^0.5*V'*w';                    % calculate and normalize the a
for i=1:3
    a(i,:)=a(i,:)/norm(a(i,:));
end
s2=a^-1*x;
if sum(s2(1,:)>0)>sum(s2(1,:)<0)
    s21=reshape(uint8(s2(1,:)),size(DW,1),size(DW,2));
else
    s21=reshape(uint8(-s2(1,:)),size(DW,1),size(DW,2));
end
figure; imshow(s21);
if sum(s2(2,:)>0)>sum(s2(2,:)<0)
    s22=reshape(uint8(s2(2,:)),size(DW,1),size(DW,2));
else
    s22=reshape(uint8(-s2(2,:)),size(DW,1),size(DW,2));
end
figure; imshow(s22);
if sum(s2(3,:)>0)>sum(s2(3,:)<0)
    s23=reshape(uint8(s2(3,:)),size(DW,1),size(DW,2));
else
    s23=reshape(uint8(-s2(3,:)),size(DW,1),size(DW,2));
end
figure; imshow(s23);
    
%%solve the original signal by FastICA Method2
W=rand(3,3);
for i=1:3
    W(:,i)=W(:,i)/norm(W(:,i));
end
W=W*(W'*W)^(-0.5);
W0=W-ones(3,3);
while sum(sum(abs(W-W0)))>0.01
    c=mean(tanh(x1'*W).*tanh(x1'*W),1)';
    W0=W;
    W=x1*tanh(x1'*W)/size(x,2)-ones(3,1)*(ones(3,1)-c)'.*W;
    for i=1:3
        W(:,i)=W(:,i)/norm(W(:,i));
    end
    W=W*(W'*W)^(-0.5);
end

A=V*D^0.5*V'*W;
for i=1:3
    A(i,:)=A(i,:)/norm(A(i,:));
end
s3=A^-1*x;
if sum(s3(1,:)>0)>sum(s3(1,:)<0)
    s31=reshape(uint8(s3(1,:)),size(DW,1),size(DW,2));
else
    s31=reshape(uint8(-s3(1,:)),size(DW,1),size(DW,2));
end
figure; imshow(s31);
if sum(s3(2,:)>0)>sum(s3(2,:)<0)
    s32=reshape(uint8(s3(2,:)),size(DW,1),size(DW,2));
else
    s32=reshape(uint8(-s3(2,:)),size(DW,1),size(DW,2));
end
figure; imshow(s32);
if sum(s3(3,:)>0)>sum(s3(3,:)<0)
    s33=reshape(uint8(s3(3,:)),size(DW,1),size(DW,2));
else
    s33=reshape(uint8(-s3(3,:)),size(DW,1),size(DW,2));
end
figure; imshow(s33);

%%we can find that the outcomings of FastICA method is better.
