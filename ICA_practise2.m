%%import the pictures, cut the original pictures

DW=imread('G:\soon\大四下\数据分析基础\DW.jpg');
DW=DW(:,:,1);
carpet=imread('G:\soon\大四下\数据分析基础\carpet.jpg');
carpet=carpet(:,:,1);
lion=imread('G:\soon\大四下\数据分析基础\lion.jpg');
lion=lion(:,:,1);

%cut out 256 16*16 images random each picture
x=zeros(1536,256);
for i=1:512
    xpos=randi(size(DW,2)-15);
    ypos=randi(size(DW,1)-15);    
    x(3*i-2,:)=reshape(DW(ypos:ypos+15,xpos:xpos+15),1,size(x,2));
    xpos=randi(size(carpet,2)-15);
    ypos=randi(size(carpet,1)-15);    
    x(3*i-1,:)=reshape(carpet(ypos:ypos+15,xpos:xpos+15),1,size(x,2));
    xpos=randi(size(lion,2)-15);
    ypos=randi(size(lion,1)-15);    
    x(3*i,:)=reshape(lion(ypos:ypos+15,xpos:xpos+15),1,size(x,2));
end

%centralize and whiten the mixed signals
x1=x-mean(x,2)*ones(1,size(x,2));
[V,D]=eig(x1*x1'/size(x,2));
x1=real(single(V*D^(-1/2)*V'*x1));

%%solve the original signal by FastICA Method2
W=rand(1536,256);
W=W*(W'*W)^(-0.5);
W=real(single(W));
W0=W-ones(1536,256);
cnum=0;
while mean(mean(abs(W-W0)))>0.01
    W0=W;
    c=mean(tanh(x1'*W).*tanh(x1'*W),1)';
    W=x1*tanh(x1'*W)/size(x,2)-ones(1536,1)*(ones(256,1)-c)'.*W;
    W=W*(W'*W)^(-0.5);
    W=real(single(W));
    cnum=cnum+1;
    if cnum>1000
        break;
    end
end

A=V*D^0.5*V'*W;
A=real(single(A));
for i=1:1536
    A(i,:)=A(i,:)/norm(A(i,:));
end
[UA,SA,VA]=svd(A);
for i=1:min(size(SA,1),size(SA,2))
    SA(i,i)=SA(i,i)^-1;
end
s=VA*SA'*UA'*x;                    %solve the independent signals
l2mod=sum(A.*A);
term=sort(l2mod);
pos=zeros(1,size(A,2));
for i=1:size(A,2)
    pos(i)=size(A,2)-find(term==l2mod(i))+1;
end
pos1=find(pos==1);
if sum(s(pos1,:)>0)>sum(s(pos1,:)<0)
    pic1=reshape(uint8(s(pos1,:)),16,16);
else
    pic1=reshape(uint8(-s(pos1,:)),16,16);
end
imshow(pic1);
pos2=find(pos==2);
if sum(s(pos2,:)>0)>sum(s(pos2,:)<0)
    pic2=reshape(uint8(s(pos2,:)),16,16);
else
    pic2=reshape(uint8(-s(pos2,:)),16,16);
end
figure;imshow(pic2);
pos3=find(pos==3);
if sum(s(pos3,:)>0)>sum(s(pos3,:)<0)
    pic3=reshape(uint8(s(pos3,:)),16,16);
else
    pic3=reshape(uint8(-s(pos3,:)),16,16);
end
figure;imshow(pic3);
