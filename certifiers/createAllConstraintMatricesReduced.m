function [As] = createAllConstraintMatricesReduced()
D=12;
e1=1; e4=2;e7=3; e2=4; e5=5; e8=6; e3=7; e6=8; e9=9;
t1=10; t2=11;t3=12;

A2 = zeros(D, D); A2(e1,e1)=1; A2(e2,e2)=1;A2(e3,e3)=1;A2(t2, t2)=-1; A2(t3, t3)=-1;
A5 = zeros(D, D); A5(e4,e4)=1; A5(e5,e5)=1;A5(e6,e6)=1;A5(t1, t1)=-1; A5(t3, t3)=-1; 
A7 = zeros(D, D); A7(e7,e7)=1; A7(e8,e8)=1;A7(e9,e9)=1;A7(t1, t1)=-1; A7(t2, t2)=-1;

A3 = zeros(D, D); A3(e1,e4)=1; A3(e2,e5)=1;A3(e3,e6)=1;A3(t1, t2)=1;               
A4 = zeros(D, D); A4(e1,e7)=1; A4(e2,e8)=1;A4(e3,e9)=1;A4(t1, t3)=1;              
A6 = zeros(D, D); A6(e4,e7)=1; A6(e5,e8)=1;A6(e6,e9)=1;A6(t2, t3)=1;             

A1 = zeros(D, D); A1(t1,t1)=1; A1(t2,t2)=1;A1(t3,t3)=1;                      

As=zeros(12, 12, 7);
As(:, :, 1)=symmetrize(A1); As(:, :, 2)=symmetrize(A2);
As(:, :, 3)=symmetrize(A3); As(:, :, 4)=symmetrize(A4);
As(:, :, 5)=symmetrize(A5); As(:, :, 6)=symmetrize(A6);
As(:, :, 7)=symmetrize(A7);

end

