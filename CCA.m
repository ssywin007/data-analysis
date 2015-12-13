%practice 50, with sample covarience matrix
x=[191 155;195 149;181 148;183 153;176 144;208 157;189 150;197 159;188 152;192 150;179 158;183 147;174 150;190 159;188 151;163 137;195 155;186 153;181 145;175 140;192 154;174 143;176 139;197 167;190 163];
y=[179 145;201 152;185 149;188 149;171 142;192 152;190 149;189 152;197 159;187 151;186 148;174 147;185 152;195 157;187 158;161 130;183 158;173 148;182 146;165 137;185 152;178 147;176 143;200 158;187 150];
s=cov([x,y]);
%block s into s11, s12, s21 and s22
s11=s(1:size(x,2),1:size(x,2));
s22=s(size(x,2)+1:size(s,2),size(x,2)+1:size(s,2));
s12=s(1:size(x,2),size(x,2)+1:size(s,2));
s21=s12';
S=s11^-1*s12*s22^-1*s21;
T=s22^-1*s21*s11^-1*s12;
[VS,DS]=eig(S);
[VT,DT]=eig(T);
%calculate the canonical variates and canonical correlation coefficient
lamda=zeros(min(size(x,2),size(y,2)),1);
a=zeros(size(x,2),min(size(x,2),size(y,2)));
b=zeros(size(y,2),min(size(x,2),size(y,2)));
[BS,PS]=sort(max(DS),'descend');
[BT,PT]=sort(max(DT),'descend');
for i=1:min(size(x,2),size(y,2))
    lamda(i)=BS(i)^0.5;
    a(:,i)=VS(:,PS(i))/(VS(:,PS(i))'*s11*VS(:,PS(i)))^0.5;
    b(:,i)=VT(:,PT(i))/(VT(:,PT(i))'*s22*VT(:,PT(i)))^0.5;
end
%U1=0.0566x1+0.0707x2;U2=-0.14x1+0.1869x2;
%V1=-0.0502y1-0.0802y2;V2=-0.1761y1+0.2621y2;
%lamda1=0.7885,lamda2=0.0537
U=(a'*(x-(ones(size(x,1),1)*mean(x)))')';
V=(b'*(y-(ones(size(y,1),1)*mean(y)))')';
UV=cov([U,V]);
%I find that the correaltion coefficient of U1 & V1 is -0.7885 and the
%correaltion coefficient of U2 & V2 is 0.0537, however, the lamdas are 0.7885 and 0.0537. It is because we solve the eigenvalue dividually.
