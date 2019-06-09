var
x1
x2
x3

varexo 
e

parameters
a11
a12
a13
a21
a22
a23
a31
a32
a33
b11
b12
b13
b21
b22
b23
b31
b32
b33
paramval
a11 = 3
a12 = 4
a13 = 5
a21 = 3
a22 = 2
a23 = 1
a31 = 6
a32 = 2
a33 = 5
b11 = 8
b12 = 7
b13 = 3
b21 = 1
b22 = 4
b23 = 2
b31 = 9
b32 = 2
b33 = 4

model
a11*x1 = b11*x1(-1) + b12*x2(-1) + 
b13*x3(-1)
 + e;
a21*x1 + a23*x3 = b23*x3(-1);
a33*x3 = b31*x1(-1) + b32*x2(-1); 
a11*x1 + a22*x2 + a33*x3 = b11*x1(-1) + b22*x2(-1) + b33*x3(-1); 

