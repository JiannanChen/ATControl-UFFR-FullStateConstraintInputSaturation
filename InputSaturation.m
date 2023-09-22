function [sys,x0,str,ts,simStateCompliance] = InputSaturation(t,x,u,flag)

switch flag,

  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;

  case 1,
    sys=mdlDerivatives(t,x,u);

  case 2,
    sys=[];

  case 3,
    sys=mdlOutputs(t,x,u);

  case 4,
    sys=[];

  case 9,
    sys=[];

  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% mdlInitializeSizes

function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 4;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 7;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

% initialize the initial conditions

x0  = [0.1, 0.2, 0.1, 0.3];

% str is always an empty matrix

str = [];

% initialize the array of sample times

ts  = [0 0];

simStateCompliance = 'UnknownSimState';


function sys=mdlDerivatives(t, x, u)
% 消防车物理参数
m = 200;
c = 0.05;
r = 0.1;
Jw = 0.005;
Jv = 10;
k = 6;
Fp = 50;
s = 1;
beta = 0.1*sin(t)*180/pi;
l = 0.3;

% 系统参数
a1 = 2*c/(m*r^2 + 2*Jw);
b1 = k*r/(m*r^2 + 2*Jw);
c1 = r*r/(m*r^2 + 2*Jw);
a2 = 2*c*l*l/(Jv*r^2 + 2*Jw*l^2);
b2 = r*l*k/(Jv*r^2 + 2*Jw*l^2);
c2 = r*r/(Jv*r^2 + 2*Jw*l^2);
g1 = 0.9;
g2 = 0.95;
B1 = b1*g1;
B2 = b2*g2;
Fdv = -Fp*cos(beta);
Fdp = -s*Fp*cos(beta)*sin(beta);

% 控制环节部分
% 非线性项
f1 = x(1)*x(1) + 2;
% 控制参数
k1 = 1000;
k2 = 100;
cv = 10;
l2 = 0.01;
gama1 = 10; 
epsilon1 = 100;
delta1 = 1;

% 限制条件
y1d = 2*sin(t);
dy1d = 2*cos(t);

F11 = 2.5;
F12 = 3;
dF11 = 0;
dF12 = 0;
F21 = 30;
F22 = 30;
dF21 = 0;
dF22 = 0;

% 变换
xi_1d = y1d * csch(F11+y1d) * csch(F12-y1d);
dxi_1d = dy1d * csch(F11+y1d) * csch(F12-y1d) ...
               + y1d * (-csch(F11+y1d)*coth(F11+y1d)*(dF11+dy1d)) * csch(F12-y1d) ...
               + y1d * csch(F11+y1d) * (-csch(F12-y1d)*coth(F12-y1d)*(dF12-dy1d));

xi_1 = x(1) * csch(F11+x(1)) * csch(F12-x(1));
xi_2 = x(2);

rho_1 = sinh(F11+x(1)) * sinh(F12-x(1));
h1 = -x(1) * (dF11*coth(F11+x(1)) + dF12*coth(F12-x(1))) / rho_1;
rho_v = ((F22-F21)*x(2) + (((F22+F21)*x(2))^2 - 2*(F22-F21)*x(2)+1)^(1/2) - 1) / (2*x(2)*x(2));
bw1 = (1 - x(1)*coth(F11+x(1)) + x(1)*coth(F12-x(1))) / rho_1;


% 误差
z1 = xi_1 - xi_1d;
z2 = xi_2 - x(3);


% 饱和输入环节
Uv = rho_v*x(2);

% 非线性项
Phi1 = bw1^2*f1^2*epsilon1 + h1^2*epsilon1 + 3*bw1^2*epsilon1 + dxi_1d^2*epsilon1;

% 辅助控制信号
alpha1 = (-1/bw1) * (k1*z1 + x(4)*z1*Phi1);
da2f = (1/l2) * (-x(3) + alpha1 / rho_v);
wv = -k2*z2 + cv*xi_2 + da2f - bw1^2*rho_v^2*z1^2*z2*epsilon1;

% 系统模型部分
sys(1) = -a1*x(1) + B1*Uv - c1*Fdv;  % 原始系统状态x1
sys(2) = -cv*x(2) + wv;  % 扩展状态
sys(3) = (1/l2) * (-x(3) + alpha1 / rho_v);  % 滤波器
sys(4) = gama1*z1^2*Phi1 - delta1*x(4); % 自适应律


function sys=mdlOutputs(t, x, u)
% 消防车物理参数
m = 200;
c = 0.05;
r = 0.1;
Jw = 0.005;
Jv = 10;
k = 6;
Fp = 50;
s = 1;
beta = 0.1*sin(t)*180/pi;
l = 0.3;

% 系统参数
a1 = 2*c/(m*r^2 + 2*Jw);
b1 = k*r/(m*r^2 + 2*Jw);
c1 = r*r/(m*r^2 + 2*Jw);
a2 = 2*c*l*l/(Jv*r^2 + 2*Jw*l^2);
b2 = r*l*k/(Jv*r^2 + 2*Jw*l^2);
c2 = r*r/(Jv*r^2 + 2*Jw*l^2);
g1 = 0.9;
g2 = 0.95;
B1 = b1*g1;
B2 = b2*g2;
Fdv = -Fp*cos(beta);
Fdp = -s*Fp*cos(beta)*sin(beta);

% 控制环节部分
% 非线性项
f1 = x(1)*x(1) + 2;
% 控制参数
k1 = 1000;
k2 = 100;
cv = 10;
l2 = 0.01;
gama1 = 10; 
epsilon1 = 100;
delta1 = 1;

% 限制条件
y1d = 2*sin(t);
dy1d = 2*cos(t);

F11 = 2.5;
F12 = 3;
dF11 = 0;
dF12 = 0;
F21 = 30;
F22 = 30;
dF21 = 0;
dF22 = 0;

% 变换
xi_1d = y1d * csch(F11+y1d) * csch(F12-y1d);
dxi_1d = dy1d * csch(F11+y1d) * csch(F12-y1d) ...
               + y1d * (-csch(F11+y1d)*coth(F11+y1d)*(dF11+dy1d)) * csch(F12-y1d) ...
               + y1d * csch(F11+y1d) * (-csch(F12-y1d)*coth(F12-y1d)*(dF12-dy1d));

xi_1 = x(1) * csch(F11+x(1)) * csch(F12-x(1));
xi_2 = x(2);

rho_1 = sinh(F11+x(1)) * sinh(F12-x(1));
h1 = -x(1) * (dF11*coth(F11+x(1)) + dF12*coth(F12-x(1))) / rho_1;
rho_v = ((F22-F21)*x(2) + (((F22+F21)*x(2))^2 - 2*(F22-F21)*x(2)+1)^(1/2) - 1) / (2*x(2)*x(2));
bw1 = (1 - x(1)*coth(F11+x(1)) + x(1)*coth(F12-x(1))) / rho_1;
% 误差
z1 = xi_1 - xi_1d;
z2 = xi_2 - x(3);
% 饱和输入环节
Uv = rho_v*x(2);
% 非线性项
Phi1 = bw1^2*f1^2*epsilon1 + h1^2*epsilon1 + 3*bw1^2*epsilon1 + dxi_1d^2*epsilon1;
% 辅助控制信号
alpha1 = (-1/bw1) * (k1*z1 + x(4)*z1*Phi1);
da2f = (1/l2) * (-x(3) + alpha1 / rho_v);
wv = -k2*z2 + cv*xi_2 + da2f - bw1^2*rho_v^2*z1^2*z2*epsilon1;

sys(1) = x(1);
sys(2) = y1d;
sys(3) = Uv;
sys(4) = -F11;
sys(5) = F12;
sys(6) = -F21;
sys(7) = F22;




