% Attempt to reproduce Fig. 8A of Levy 2018

% 1. particle and surface properties
E_p   = 700e9; % Pa
Y_p   = 1500.4e9;
rho_p = 3.89e3; % kg/m3
nv_p  = 0.29;
H_p   = 2550e9; % MPa

E_w   = 200e9;  % Pa
Y_w   = 205e6;
rho_w = 8.03e3; % kg/m3
nv_w  = 0.3;
H_w   = 2.4e9; % Pa 

% 2. derived properties
E_star = (1-nv_p^2)/E_p + (1-nv_w^2)/E_w;
E_star = E_star^-1;

v_yw   = 5.052*(Y_w^5/E_star^4/rho_w)^0.5;
v_iw   = 0.02*v_yw*(E_star/Y_w)^2;

% model parameters
Ch = 66;

% e_n results
alpha = 60/180*pi;
N   = 200;
v_i = linspace(0,350,N)';
e_n = v_i;

for i=1:N
  vi_vyw_ch = v_i(i) / v_yw*sin(alpha)^3.5/Ch;
  if(vi_vyw_ch<100)
    e_n(i) = (vi_vyw_ch)^-0.091;
  elseif (vi_vyw_ch < v_iw/v_yw)
    e_n(i) = 2.08*(vi_vyw_ch)^-0.25;
  else
    e_n(i) = 0.78*(vi_vyw_ch/(E_star/Y_w))^-0.5;
  end
end

levy = readmatrix('Input_files/levy_60_degrees.dat');
plot(v_i,e_n,'m',LineWidth=2); hold on
plot(levy(:,1),levy(:,2),'ro',LineWidth=2);
legend('Eq. 21, E_p = 700 GPa, \nu = 0.29','Levy - Fig. 8a');
xlabel('impact velocity, v_i (m/s)');
ylabel('normal coefficient of restitution, e_n');
title('hardened steel ball impacting obliquely on SS 301 (Fig. 8A, Levy 2018)');
grid on
axis([0 350 0.0 1.0]);
% text(200,0.8,'impact angle = 15 degree');
hold on

