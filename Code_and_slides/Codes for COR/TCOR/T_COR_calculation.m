clear
clc
%% Open pair property file
fileID1 = fopen("Input_files/file1.txt","r");
formatSpec1 = ['%s' '%f' '%f' '%f'];
pair_property = textscan(fileID1,formatSpec1,'headerlines',1);
fclose(fileID1);
%% .........User gives 1st input: (wall,particle) pair...........%
pair = 'gamma-Al2O3,Steel';
%% TCOR- Code extracts required material properties from database
% gamma-Al2O3, Steel
for i = 1:length(pair_property{1,1})
    if pair == string(pair_property{1,1}(i))
        c3 = pair_property{1,3}(i);
        c4 = pair_property{1,4}(i);
        T1 = readtable('Input_files/TCOR_impact_angle_gamma_Al2O3.csv');
        x1.alpha = T1.Var1;
        x1.nu_w = 0.29;
        x1.mu = 0.131;
        x1.ncor_expt = 0.735;
        y1 = T1.Var2; % Experimental data of TCOR
        p1 = [1,-1]; % Initial guess for c3 = p(1), c4 = p(2)
        [bestP1,bestErr1] = best_fit('TCOR_Err',p1,x1,y1);
        err_fitting_1 = [abs(bestP1(1)-c3)/abs(c3),abs(bestP1(2)-c4)/abs(c4)]*100;
        my_tcor = t_cor_wall(x1,bestP1);
        levy_tcor = t_cor_wall(x1,[c3,c4]);
        error_my_tcor = (my_tcor-y1')./y1'*100;
        error_levy_tcor = (levy_tcor-y1')./y1'*100;
        figure(1)
        subplot(2,1,1)
        plot(x1.alpha,y1,x1.alpha,my_tcor,x1.alpha,levy_tcor,'o',LineWidth=2);
        xlabel('impact angle');
        ylabel('TCOR');
        legend('Expt.','My TCOR','Levy TCOR');
        title('\gamma-Al2O3,Steel');
        subplot(2,1,2)
        plot(x1.alpha,error_my_tcor,x1.alpha,error_levy_tcor,LineWidth=2);
        xlabel('impact angle');
        ylabel('Deviation from expt.value (%)');
        legend('My TCOR','Levy TCOR');
        title('\gamma-Al2O3,Steel');
        break
    end
end
% Sodium-benzoate, Steel
pair = 'Sodium-benzoate,Steel';
for i = 1:length(pair_property{1,1})
    if pair == string(pair_property{1,1}(i))
        c3 = pair_property{1,3}(i);
        c4 = pair_property{1,4}(i);
        % sodium-benzoate on steel
        T2 = readtable('Input_files/TCOR_impact_angle_sodium_benzoate.csv');
        x2.alpha = T2.Var1;
        x2.nu_w = 0.29;
        x2.mu = 0.153;
        x2.ncor_expt = 0.532;
        y2 = T2.Var2; % Experimental data of TCOR
        p2 = [1,-1]; % c3 = p(1), c4 = p(2)
        [bestP2,bestErr2] = best_fit('TCOR_Err',p2,x2,y2);
        err_fitting_2 = [abs(bestP2(1)-c3)/abs(c3),abs(bestP2(2)-c4)/abs(c4)]*100;
        my_tcor = t_cor_wall(x2,bestP2);
        levy_tcor = t_cor_wall(x2,[c3,c4]);
        error_my_tcor = (my_tcor-y2')./y2'*100;
        error_levy_tcor = (levy_tcor-y2')./y2'*100;
        figure(2)
        subplot(1,2,1)
        plot(x2.alpha,y2,x2.alpha,my_tcor,x2.alpha,levy_tcor,'o',LineWidth=2);
        xlabel('impact angle');
        ylabel('TCOR');
        legend('Expt.','My TCOR','Levy TCOR');
        title('Sodium-benzoate,Steel');
        subplot(1,2,2)
        plot(x2.alpha,error_my_tcor,x2.alpha,error_levy_tcor,LineWidth=2);
        xlabel('impact angle');
        ylabel('Deviation from expt.value');
        legend('My TCOR','Levy TCOR');
        title('Sodium-benzoate,Steel');
        break
    end
end
% Zeolite-4A, Steel 
pair = 'Zeolite_4A,Steel';
for i = 1:length(pair_property{1,1})
    if pair == string(pair_property{1,1}(i))
        c3 = pair_property{1,3}(i);
        c4 = pair_property{1,4}(i);
        % Zeolite-4A on steel
        T3 = readtable('Input_files/TCOR_impact_angle_zeolite_4A.csv');
        x3.alpha = T3.Var1;
        x3.nu_w = 0.29;
        x3.mu = 0.114;
        x3.ncor_expt = 0.653;
        y3 = T3.Var2; 
        p3 = [1,-1]; % c3 = p(1), c4 = p(2)
        [bestP3,bestErr3] = best_fit('TCOR_Err',p3,x3,y3);
        err_fitting_3 = [abs(bestP3(1)-c3)/abs(c3),abs(bestP3(2)-c4)/abs(c4)]*100;
        my_tcor = t_cor_wall(x3,bestP3);
        levy_tcor = t_cor_wall(x3,[c3,c4]);
        error_my_tcor = (my_tcor-y3')./y3'*100;
        error_levy_tcor = (levy_tcor-y3')./y3'*100;
        figure(3)
        subplot(1,2,1)
        plot(x3.alpha,y3,x3.alpha,my_tcor,x3.alpha,levy_tcor,'o',LineWidth=2);
        xlabel('impact angle');
        ylabel('TCOR');
        legend('Expt.','My TCOR','Levy TCOR');
        title('Zeolite 4A,Steel');
        subplot(1,2,2)
        plot(x3.alpha,error_my_tcor,x3.alpha,error_levy_tcor,LineWidth=2);
        xlabel('impact angle');
        ylabel('Deviation from expt.value');
        legend('My TCOR','Levy TCOR');
        title('Zeolite 4A,Steel');
        break
    end
end
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
%% Error function for TCOR
function err= TCOR_Err(p,x,y)
    pred = t_cor_wall(x,p);
    err = sum((pred(:)-y(:)).^2);
end
%% Wall tangential coefficient of restitution
%......Impact conditions..........%
% 1. impact angle, alpha
% 2. normal impact velocity, v_ni
%......Material properties..........%
% 3. Coefficient for energy loss as heat, C_h
% 4. yield stress of wall, Y_w
% 5. Poisson's ratio of wall, nu_w
% 6. Poisson's ratio of particle, nu_p
% 7. Wall's Young's modulus, E_p
% 8. Particle's Young's modulus, E_w
% 9. Wall's material density, rho_w
% 10. Sliding friction coefficient, mu
%.....Other parameters.........
% 11. Coefficients c3,c4 
% 12. Experimental wall N_COR, n_cor_expt
function T_COR_wall = t_cor_wall(x,p)
    k = 2*(1-x.nu_w)/(2-x.nu_w);
    theta_cr = (7*k-1)/k; % critical dimensionless impact angle
    e_n = x.ncor_expt;
    theta = 2*tan(deg2rad(x.alpha))'/(x.mu*(1+e_n)); % dimensionless impact angle
    fac1 = tanh(p(1)+p(2)*theta_cr); 
    fac2 = tanh(p(1)+p(2).*theta);
    fac3 = tanh(p(1));
    impl_ratio = zeros(1,1);
    for i = 1:length(theta)
        if theta(i) < theta_cr
            impl_ratio(i) = x.mu.*(1+(-fac1+fac2(i))/(fac1-fac3));
        else
            impl_ratio(i) = x.mu;
        end
    end
    T_COR_wall = 1 - (2./theta).*(impl_ratio./x.mu);    
end