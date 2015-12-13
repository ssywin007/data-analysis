x11=[0.71 1.66 2.01 2.16 2.42 2.42 2.56 2.6 3.31 3.64 3.74 3.74 4.39 4.5 5.07 5.26 8.15 8.24];
x12=[2.2 2.93 3.08 3.49 4.11 4.95 5.16 5.54 5.68 6.25 7.25 7.9 8.85 11.96 15.54 15.89 18.3 18.59];
x13=[2.25 3.93 5.08 5.82 5.84 6.89 8.5 8.56 9.44 10.52 13.46 13.57 14.76 16.41 16.96 17.56 22.82 29.13];
x21=[2.2 2.69 3.54 3.75 3.83 4.08 4.27 4.53 5.32 6.18 6.22 6.33 6.97 6.97 7.52 8.36 11.65 12.45];
x22=[4.04 4.16 4.42 4.93 5.49 5.77 5.86 6.28 6.97 7.06 7.78 9.23 9.34 9.91 13.46 18.4 23.89 26.39];
x23=[2.71 5.43 6.38 6.38 8.32 9.04 9.56 10.01 10.08 10.62 13.8 15.99 17.9 18.25 19.32 19.87 21.6 22.25];
xmean=mean([mean(x11) mean(x12) mean(x13) mean(x21) mean(x22) mean(x23)]);
%the division of varience
SST=(x11-xmean)*(x11-xmean)'+(x12-xmean)*(x12-xmean)'+(x13-xmean)*(x13-xmean)'+(x21-xmean)*(x21-xmean)'+(x22-xmean)*(x22-xmean)'+(x23-xmean)*(x23-xmean)';
SSE=(x11-mean(x11))*(x11-mean(x11))'+(x12-mean(x12))*(x12-mean(x12))'+(x13-mean(x13))*(x13-mean(x13))'+(x21-mean(x21))*(x21-mean(x21))'+(x22-mean(x22))*(x22-mean(x22))'+(x23-mean(x23))*(x23-mean(x23))';
%factor A is the kind of ferric ion
Amean1=mean([mean(x11) mean(x12) mean(x13)]);
Amean2=mean([mean(x21) mean(x22) mean(x23)]);
SSA=(Amean1-xmean)^2*3*size(x11,2)+(Amean2-xmean)^2*3*size(x11,2);
%factor B is the dosage
Bmean1=mean([mean(x11) mean(x21)]);
Bmean2=mean([mean(x12) mean(x22)]);
Bmean3=mean([mean(x13) mean(x23)]);
SSB=(Bmean1-xmean)^2*2*size(x11,2)+(Bmean2-xmean)^2*2*size(x11,2)+(Bmean3-xmean)^2*2*size(x11,2);
SSAB=size(x11,2)*(mean(x11)-Amean1-Bmean1+xmean)^2+size(x12,2)*(mean(x12)-Amean1-Bmean2+xmean)^2+size(x13,2)*(mean(x13)-Amean1-Bmean3+xmean)^2+size(x21,2)*(mean(x21)-Amean2-Bmean1+xmean)^2+size(x22,2)*(mean(x22)-Amean2-Bmean2+xmean)^2+size(x23,2)*(mean(x23)-Amean2-Bmean3+xmean)^2;
if SST-SSA-SSB-SSAB-SSE>0.01
    error('THE DIVISION OF VARIENCE HAS A MISTAKE');
end
MSE=SSE/(2*3*(size(x11,2)-1));
MSA=SSA/(2-1);
MSB=SSB/(3-1);
MSAB=SSAB/(2-1)/(3-1);
%Ftest
FA=MSA/MSE;
FB=MSB/MSE;
FAB=MSAB/MSE;
pA=1-fcdf(FA,2-1,2*3*(size(x11,2)-1));
pB=1-fcdf(FB,3-1,2*3*(size(x11,2)-1));
pAB=1-fcdf(FAB,(2-1)*(3-1),2*3*(size(x11,2)-1));
%we care the interaction impact most
if pAB<0.05
    fprintf('WHEN ALPHA=0.05, THE INTERACTION BETWEEN THE KIND OF FERRIC ION AND THE DOSAGE HAS A NOTABLE IMPACT ON THE RESIDUAL DOSAGE.')
else
    fprintf('WHEN ALPHA=0.05, THE INTERACTION BETWEEN THE KIND OF FERRIC ION AND THE DOSAGE HAS NO NOTABLE IMPACT ON THE RESIDUAL DOSAGE.')
end
if pA<0.05
    fprintf('WHEN ALPHA=0.05, THE KIND OF FERRIC ION HAS A NOTABLE IMPACT ON THE RESIDUAL DOSAGE.')
else
    fprintf('WHEN ALPHA=0.05, THE KIND OF FERRIC ION HAS NO NOTABLE IMPACT ON THE RESIDUAL DOSAGE.')
end
if pB<0.05
    fprintf('WHEN ALPHA=0.05, THE DOSAGE HAS A NOTABLE IMPACT ON THE RESIDUAL DOSAGE.')
else
    fprintf('WHEN ALPHA=0.05, THE DOSAGE HAS NO NOTABLE IMPACT ON THE RESIDUAL DOSAGE.')
end
%pp=anova2([x11',x21';x12',x22';x13',x23'],size(x11,2)); to ensure our calculation is correct.
