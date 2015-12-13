%read umist_cropped first
p=size(facedat{1,1},3);
n=size(facedat{1,1},1)*size(facedat{1,1},2);
face=zeros(p,n);
for i=1:p
    face(i,:)=reshape(facedat{1,1}(:,:,i),n,1);
end
face=single(face);
meanface=mean(face);
stdface=std(face);
%standardlize
for i=1:p
    face(i,:)=(face(i,:)-meanface)./stdface;
end
c=cov(face');
[v,d]=eig(c);
%6 main factors
mainf=zeros(p,6);
for i=1:6
    mainf(:,i)=v(:,p+1-i)/(v(:,p+1-i)'*v(:,p+1-i));
end
mainface=mainf'*face;
for i=1:6
    mainface(i,:)=mainface(i,:).*stdface+meanface;
end
mainface=uint8(mainface);
for i=1:6
    term=reshape(mainface(i,:),size(facedat{1,1},1),size(facedat{1,1},2));
    figure;imshow(term);
end
meanface=reshape(meanface,size(facedat{1,1},1),size(facedat{1,1},2));
meanface=uint8(meanface);
figure;imshow(meanface);
