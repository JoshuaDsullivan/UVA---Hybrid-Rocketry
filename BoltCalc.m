%Bolt Calc Imperial
clc; clear; close all
%% Input
n=[4 6 8 10 12 13 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50]; %# of bolt on one end
d=[1/8 3/16 1/4]; %in %diameter of bolt
dis=.3; %in %bolt hole distance from edge of casing
MinFS=[2 3 4]; %Minimum Factor of Safety
%% Setup
P=809; %psi
d_cc=4.75%in
t=.4; %in
t_c=.125;
sigma_y_c=30000; %psi %Aluminum 6061T6
sigma_y_t=45000; %psi %Aluminum 6061T6
tau=0.577*sigma_y_t; %psi %Conservative Conversion
sigma_y_t_st=170000; %psi %Alloy Steel, update from supplier
tau_st=0.577*sigma_y_t_st; %psi %Conservative Conversion
sigma_b=[];
sigma_t=[];
sigma_s=[];
FS_1=zeros(length(n),length(d));
FS_2=zeros(length(n),length(d));
FS_3=zeros(length(n),length(d));
%% Results
for i=1:length(n)
    for j=1:length(d)
        sigma_b=(P*pi*(d_cc/2)^2)/(n(i)*d(j)*t);
        sigma_t=(P*pi*(d_cc/2)^2)/(2*n(i)*dis*t);
        sigma_s=(P*pi*(d_cc/2)^2)/(n(i)*pi*((d(j)/2)^2));
        FS_1(i,j)=sigma_y_c/sigma_b;
        FS_2(i,j)=tau/sigma_t;
        FS_3(i,j)=tau_st/sigma_s;
    end
end
%% Plot
%% Plot
[n_plot,d_plot]=ndgrid(n,d);
figure;
scatter3(n_plot,d_plot,FS_1,50,'filled')
hold on
scatter3(n_plot,d_plot,FS_2,50,'filled')
scatter3(n_plot,d_plot,FS_3,50,'filled')

Z = 2 * ones(size(n_plot)); % Factor of Safety Plane
surf(n_plot, d_plot, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'cyan')

xlabel("# of Bolts")
ylabel("Bolt Diameter (in)")
zlabel("Factor of Safety")
title("Factor of Safety for 3 Failure Modes")
hold off

iter=0;
for k=1:length(MinFS)
    for l=1:length(n)
        for m=1:length(d)
            if FS_1(l,m)>=MinFS(k) && FS_2(l,m)>=MinFS(k) && FS_3(l,m)>=MinFS(k)
                iter = iter+1;
                Result(iter,:)=[MinFS(k),n(l),d(m), FS_1(l,m), FS_2(l,m), FS_3(l,m)];
            end
        end
    end
end
Result=array2table(Result);
Result.Properties.VariableNames = ["Minimum Factor of Safety","Number of Bolts","Bolt Diameter (in)", "Bearing", "Tear", "Shear"];
disp("Usable Options:")
disp(Result)
