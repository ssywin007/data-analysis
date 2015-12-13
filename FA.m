%main factor method
X=[77 82 67 67 81;75 73 71 66 81;63 63 65 70 63;51 67 65 65 68;62 60 58 62 70;52 64 60 63 54;50 50 64 55 63;31 55 60 57 73;44 69 53 53 53;62 46 61 57 45;44 61 52 62 46;12 58 61 63 67;54 49 56 47 53;44 56 55 61 36;46 52 65 50 35;30 69 50 52 45;40 27 54 61 61;36 59 51 45 51;46 56 57 49 32;42 60 54 49 33;23 55 59 53 44;41 63 49 46 34; 63 78 80 70 81;55 72 63 70 68;53 61 72 64 73;59 70 68 62 56;64 72 60 62 45;55 67 59 62 44;65 63 58 56 37;60 64 56 54 40;42 69 61 55 45;31 49 62 63 62;49 41 61 49 64;49 53 49 62 47;54 53 46 59 44;18 44 50 57 81;32 45 49 57 64;46 49 53 59 37;31 42 48 54 68;56 40 56 54 35;45 42 55 56 40;40 63 53 54 25;48 48 49 51 37;46 52 53 41 40];
%standardize the X
Xnorm=(X-ones(size(X,1),1)*mean(X))./(ones(size(X,1),1)*std(X));
s=cov(Xnorm);
[v,d]=eig(s);
%find out the main factors
m=0;
p0=0;
while p0<0.8
    p0=p0+d(size(d,1)-m,size(d,1)-m)/sum(sum(d));
    m=m+1;
end
%m public factors
A=zeros(size(d,1),m);
for i=1:m
    A(:,i)=d(size(d,1)+1-i,size(d,1)+1-i)^0.5*v(:,size(d,1)+1-i);
end
%calculate the varience  of A
VA=V(A);
%orthogonal rotation to maximize the varience
if m==1
    fprintf('no need to maximize the var.')
elseif m==2
        [A(:,1),A(:,2)]=maxvar(A(:,1),A(:,2));
else
    VA0=VA-1;
    while abs(VA-VA0)>0.0001
        for i=1:m
            for j=i+1:m
                [A(:,i),A(:,j)]=maxvar(A(:,i),A(:,j));
            end
        end
        VA0=VA;
        VA=V(A);
    end
end
%calculate the cimmunalities
h=sum(A.*A,2);
%special factors
D=eye(5)-diag(h);
%calculate the varience contribution of public factors
q=sum(A.*A,1);
B=X*s^-1;
%calculate the score of factors
F=A'*s^-1*Xnorm';

%A suggests that x1 is mainly attributed to F3, x2 is mainly attributed to F2,  x3 is mainly attributed to F1 and F3, x4 and x5 are mainly attributed to F1.
%h suggests that x1, x2 and x5 are more relied on the public factors than x3 and x4.
%q suggests F1 contribute more to the varience of X
