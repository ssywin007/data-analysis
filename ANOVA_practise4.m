A=[35.6 0.88 17.5 3.5 6.15 5.22 1.81;34.1 0.88 18.0 3.7 6.17 5.36 1.35;32.65 0.88 16.5 3.5 6.11 5.26 1.6];
B=[39.7 0.88 21 6.8 5.82 3.4 9.1; 40.7 0.88 22 6.2 5.86 3.46 8.59;42.55 0.88 21 6.1 5.88 3.43 9.28];
C=[35.3 1.31 14 2.3 4.55 5.21 9.33;34.6 1.31 14 1.7 4.6 5.19 9.21;36.15 1.31 13.5 2.3 4.59 4.96 9.45];
D=[37.8 0.44 20 5.8 4.17 6.68 5.09;40.55 0.44 19 6.4 4.21 6.64 5.22;41.3 0.44 19.5 6.2 4.16 6.63 5.21];
E=[35.85 0.44 17 4.9 3.79 5.34 7.85;36.85 0.44 17 5 3.82 5.3 8.18;36.45 0.44 18.5 4.8 3.82 5.38 7.74];
meanA=mean(A);
meanB=mean(B);
meanC=mean(C);
meanD=mean(D);
meanE=mean(E);
meanx=mean([A;B;C;D;E]);
%intra-group deviation matrix
varA1=(A-ones(size(A,1),1)*meanA)'*(A-ones(size(A,1),1)*meanA);
varB1=(B-ones(size(B,1),1)*meanB)'*(B-ones(size(B,1),1)*meanB);
varC1=(C-ones(size(C,1),1)*meanC)'*(C-ones(size(C,1),1)*meanC);
varD1=(D-ones(size(D,1),1)*meanD)'*(D-ones(size(D,1),1)*meanD);
varE1=(E-ones(size(E,1),1)*meanE)'*(E-ones(size(E,1),1)*meanE);
%inter-group deviation matrix
varA2=(mean(A)-meanx)'*(mean(A)-meanx);
varB2=(mean(B)-meanx)'*(mean(B)-meanx);
varC2=(mean(C)-meanx)'*(mean(C)-meanx);
varD2=(mean(D)-meanx)'*(mean(D)-meanx);
varE2=(mean(E)-meanx)'*(mean(E)-meanx);
varT=([A;B;C;D;E]-ones(size([A;B;C;D;E],1),1)*meanx)'*([A;B;C;D;E]-ones(size([A;B;C;D;E],1),1)*meanx);
varA=varA1+varB1+varC1+varD1+varE1;
varB=size(A,1)*(varA2+varB2+varC2+varD2+varE2);
LAMDA=det(varA)/det(varA+varB);
%LAMDA=0.  Because the second feature is always the same in each group.
