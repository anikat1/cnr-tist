syms v1 v2 real;
V = [0 v1;v2 0];
x11 = 0.4423294049321163;
x12 = 0.47224497178589264;
x21 = 0.4969854274748294;
x22 = 0.5427998723063714;
x31 = 0.9467733989027951;
x32 = 0.9532423754168234;
v1 = 0.13596;
v2 = 0.14079;
X = [x11 x12;x21 x22;x31 x32];
Y1 = (x11-x12*v2)^2+(x21-x22*v2)^2+(x31-x32*v2)^2+(x12-x11*v1)^2+(x22-x21*v1)^2+(x32-x31*v1)^2;
Y2 = (x12^2*v2-x11*x12+x22^2*v2-x21*x22+x32^2*v2-x31*x32);
Y3 = (x11^2*v1-x11*x12+x21^2*v1-x21*x22+x31^2*v1-x32*v2);

dfdv1 = -0.5*Y1^(-3/2)*Y3*(x11^2+x21^2+x31^2);
dfdv2 = -0.5*Y1^(-3/2)*Y2*(x12^2+x22^2+x32^2);
dfdv12 = -0.5*Y1^(-3/2)*Y2*Y3;
H = [dfdv1 dfdv12;dfdv12 dfdv2]
%f = .5 * norm(X' - X'*V,'fro')^2 + -1.00*norm(V,1) + 2.00*l2l1norm(VR);
%H = hessian(f,[v12,v13,v21,v23,v31,v32]);
%size(H)

eigen = eig(vpa(H))
d= det(H)
%ans = min(eigen)

