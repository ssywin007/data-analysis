%run umist_cropped.m to load 1st people's umist_cropped photos
%p-number of factors-=38
p=size(facedat{1,1},3);
%n-number of observations=10304
n=size(facedat{1,1},1)*size(facedat{1,1},2);
face=zeros(p,n);
%strentch the photo matrix to a vector
for i=1:p
    face(i,:)=reshape(facedat{1,1}(:,:,i),n,1);
end
face=single(face);
meanface=(mean(face));
%substract the mean of 38 photos
for i=1:p
    face(i,:)=face(i,:)-meanface;
end
c=cov(face');
[v,d]=eig(c);
%6 main factors
mainf=zeros(p,6);
for i=1:6
    mainf(:,i)=v(:,p+1-i)/(v(:,p+1-i)'*v(:,p+1-i));
end
%calculate the main factor face
mainface=mainf'*face;
for i=1:6
    mainface(i,:)=mainface(i,:)+meanface;
end
mainface=uint8(mainface);
%represent 6 main factor faces
for i=1:6
    term=reshape(mainface(i,:),size(facedat{1,1},1),size(facedat{1,1},2));
    figure;imshow(term);
end
%represent the mean faces
meanface=reshape(meanface,size(facedat{1,1},1),size(facedat{1,1},2));
meanface=uint8(meanface);
figure;imshow(meanface);
