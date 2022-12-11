clear
clc
%% Open pair property file
fileID1 = fopen("Input_files/file1.txt","r");
formatSpec1 = ['%s' '%f' '%f' '%f'];
pair_property = textscan(fileID1,formatSpec1,'headerlines',1);
fclose(fileID1);
%% .........User gives 1st input: (wall,particle) pair...........%
pair = 'Alumina,Aluminium';
%% Wall NCOR- Code extracts required material properties from database
for i = 1:length(pair_property{1,1})
    if pair == string(pair_property{1,1}(i))
        dummy_string = string(pair_property{1,1}(i));
        % Alumina on Aluminium alloy
        T1 = readtable('Input_files/NCOR_vs_impact_vel_normal_impact_alumina_on_Al.csv');
        x1.alpha = 90;
        x1.v_i = T1.Var1;
        y1 = T1.Var2; % Experimental data of NCOR
        break
    end
end
dummy_string = split(dummy_string,',');
wall_string = dummy_string(2);
particle_string = dummy_string(1);
% Open material property file
fileID2 = fopen("Input_files/file2.txt","r");
formatSpec2 = ['%s' '%f' '%f' '%f' '%f' '%f' '%f' '%f'];
material_property = textscan(fileID2,formatSpec2,'headerlines',1);
fclose(fileID2);
for i = 1:length(material_property{1,1})
    if wall_string == string(material_property{1,1}(i))
        x1.Y_w = material_property{1,2}(i)*1e6;
        x1.nu_w = material_property{1,3}(i);
        x1.E_w = material_property{1,4}(i)*1e9;
        x1.rho_w = material_property{1,5}(i);
        Ch_1 = material_property{1,7}(i);
        break
    end
end
for i = 1:length(material_property{1,1})
    if particle_string == string(material_property{1,1}(i))
        x1.nu_p = material_property{1,3}(i);
        x1.E_p = material_property{1,4}(i)*1e9;
        break
    end
end
C_h = 1; % Initial guess for C_h
[bestC_h,bestErr1] = best_fit('NCOR_Err',C_h,x1,y1);
err_fitting_1 = (abs(bestC_h-Ch_1)/abs(Ch_1))*100;
my_ncor = n_cor_wall(x1,bestC_h);
levy_ncor = n_cor_wall(x1,Ch_1);
error_my_ncor = (my_ncor-y1')./y1'*100;
error_levy_ncor = (levy_ncor-y1')./y1'*100;
figure(1)
subplot(2,1,1)
plot(x1.v_i,y1,x1.v_i,my_ncor,x1.v_i,levy_ncor,'o',LineWidth=2);
xlabel('impact velocity');
ylabel('NCOR');
legend('Expt.','My NCOR','Levy NCOR');
title('Alumina,Aluminium');
subplot(2,1,2)
plot(x1.v_i,error_my_ncor,x1.v_i,error_levy_ncor,LineWidth=2);
xlabel('impact velocity');
ylabel('Deviation from expt.value (%)');
legend('My NCOR','Levy NCOR');
title('Alumina,Aluminium');

%% Optimization function
%......Input arguments.......%
% funName - name of function to be optimized, as a string
% p - vector of parameters to be found
% x - available input parameters to TCOR or NCOR correlations
% y - Expt. data of T-COR or N-COR
function [params,new_err] = best_fit(funName,p,x,y)
    optimize_f = str2func(funName);
    err = @(params) optimize_f(params,x,y);
    options = optimset('MaxFunEvals', 1E7, 'MaxIter', 1E7);
    params = fminsearch(err,p,options);
    new_err = optimize_f(params,x,y);
end
%% Error function for NCOR
function err_ncor = NCOR_Err(C_h,x,y)
    pred = n_cor_wall(x,C_h);
    err_ncor = sum((pred(:)-y(:)).^2);
end
%% Wall normal coefficient of restitution
%......Impact conditions..........%
% 1. impact angle, alpha
% 2. normal impact velocity, v_ni
%......Material properties..........%
% 3. Coefficient for energy loss as heat, C_h
% 4. yield stress of wall, Y_w
% 5. Poisson's ratio of wall, nu_w
% 6. Poisson's ratio of particle, nu_p
% 7. Wall's Young's modulus, E_w
% 8. Particle's Young's modulus, E_p
% 9. Wall's material density, rho_w
function N_COR_wall = n_cor_wall(x,C_h)
    E = (1-x.nu_p^2)/x.E_p + (1-x.nu_w^2)/x.E_w; 
    E_star = 1/E; % representative Young's modulus
    v_y = 5.5052*(x.Y_w^5/(E_star^4*x.rho_w))^0.5;
    v_iw_star = 0.02*v_y*(E_star/x.Y_w)^2;
    fac1 = (sin(deg2rad(x.alpha))).^3.5*x.v_i/(C_h*v_y);
    fac2 = v_iw_star/(C_h*v_y);
    for i = 1:length(x.v_i)
        if fac1(i) <= 100
            N_COR_wall(i) = fac1(i)^-0.091;
        elseif fac2 >= fac1(i) && fac1(i) > 100
            N_COR_wall(i) = 2.08*fac1(i)^-0.25;
        elseif fac1(i) > fac2
            N_COR_wall(i) = 0.78*(fac1(i)/(E_star/x.Y_w))^-0.5;
%         else
%             N_COR_wall = 1; % dummy value for dummy inputs
        end
    end
end