y1=[7.6 8.2 6.8 5.8 6.9 6.6 6.3 7.7 6];
y2=[6.7 8.1 9.4 8.6 7.8 7.7 8.9 7.9 8.3 8.7 7.1 8.4];
y3=[8.5 9.7 10.1 7.8 9.6 9.5];
ymean=(sum(y1)+sum(y2)+sum(y3))/(size(y1,2)+size(y2,2)+size(y3,2));
%the division of varience
SST=(y1-ymean)*(y1-ymean)'+(y2-ymean)*(y2-ymean)'+(y3-ymean)*(y3-ymean)';
SSE=(y1-mean(y1))*(y1-mean(y1))'+(y2-mean(y2))*(y2-mean(y2))'+(y3-mean(y3))*(y3-mean(y3))';
SSA=(mean(y1)-ymean)^2*size(y1,2)+(mean(y2)-ymean)^2*size(y2,2)+(mean(y3)-ymean)^2*size(y3,2);
if SST-SSA-SSE>0.01
    error('THE DIVISION OF VARIENCE HAS A MISTAKE');
end
%the degrees of freedom are a-1 for SSA, n-a for SSE.
MSA=SSA/(3-1);
MSE=SSE/(size(y1,2)+size(y2,2)+size(y3,2)-3);
%Ftest
F=MSA/MSE;
p=1-fcdf(F,3-1,size(y1,2)+size(y2,2)+size(y3,2)-3);
if p<0.05
    fprintf('WHEN ALPHA=0.05, THE INVESTMENT OF RESEARCH HAS A NOTABLE IMPACT ON THE ENHANCEMENT OF MANUFACTURE ABILITY.')
else
    fprintf('WHEN ALPHA=0.05, THE INVESTMENT OF RESEARCH HAS NO NOTABLE IMPACT ON THE ENHANCEMENT OF MANUFACTURE ABILITY.')
end
%pp=anova1([y1,y2,y3],[ones(1,size(y1,2)),2*ones(1,size(y2,2)),3*ones(1,size(y3,2))]);to ensrue our calculation is correct.
