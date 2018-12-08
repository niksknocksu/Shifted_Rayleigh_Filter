function func_NewNPOL_SRF(sheet,MC)

% Use excel file NPOL_new.xlsx to generate observer data and target data
file_in = 'NPOL_new.xlsx';
%sheet = 12 ;
n = 4;  % Dimension on state vector
A = xlsread(file_in,sheet,'');
T =1;
%kmax = numel(A.data.Sheet1(:,1));
x_obs(1,:) = A(:,1)';
x_obs(2,:) = A(:,2)';
x_obs(3,:) = A(:,3)';
x_obs(4,:) = A(:,4)';

x_tar(1,:) = A(:,5)';
x_tar(2,:) = A(:,6)';
x_tar(3,:) = A(:,7)';
x_tar(4,:) = A(:,8)';

F =[1 0 T 0;...
    0 1 0 T;...
    0 0 1 0;...
    0 0 0 1]; % State Transition matrix
q=[T^3/3 0 T^2/2 0;...
   0 T^3/3 0 T^2/2;...
   T^2/2 0 T 0;...
   0 T^2/2 0 T];
Q=((1.944E-6)/216000)*q; %Process noise covariance
H = [1,0,0,0;...
     0,1,0,0];
%{
[theta_obs,rho_obs] = cart2pol(x_obs(1,:),x_obs(2,:));
[theta_tar,rho_tar] = cart2pol(x_tar(1,:),x_tar(2,:));
polarplot(theta_obs,rho_obs,'r');
hold on
polarplot(theta_tar,rho_tar,'b');
%}
kmax = numel(x_obs(1,:));
sigma_theta = 2*pi/180;
R = sigma_theta^2;
%MC = 1000;
track_loss = 0;
for mm = 1:MC
    fprintf('Run %d/%d: ' ,mm,MC);
    % Relative state vector
    x_rel(:,1) = x_tar(:,1)-x_obs(:,1);
    
   for k = 2:kmax
    x_rel(:,k) = (F*x_rel(:,k-1)) + mvnrnd([0,0,0,0],Q)' - [x_obs(1,k) - x_obs(1,k-1) - (T*x_obs(3,k-1));...
                                                x_obs(2,k) - x_obs(2,k-1) - (T*x_obs(4,k-1));...
                                                x_obs(3,k) - x_obs(3,k-1);...
                                                x_obs(4,k) - x_obs(4,k-1)];
   end
   x_tar = x_rel + x_obs;
    init_r = norm([x_rel(1,1),x_rel(2,1)]);
    tar_speed = norm([x_tar(3,1),x_tar(4,1)]);
    
    % Measurements
    z = zeros(1,kmax);
    for i = 1:kmax
        z(i) = atan2(x_rel(1,i),x_rel(2,i)) + mvnrnd(0,R);
    end
    
    % Initialization 
    x_rel_obs = zeros(4,kmax);
    sigma_r = 1;
    init_r = init_r + sigma_r*randn(1);
    init_sp = tar_speed;
    sigma_sp = 1*0.0308667/60; %kms/sec
    init_sp = init_sp + sigma_sp*randn(1);
    init_co = z(1);
    sigma_c = pi/sqrt(12);
    x_rel_obs(:,1) = [init_r*sin(z(1));...
                      init_r*cos(z(1));...
                      (init_sp*sin(init_co))-x_obs(3,1);...
                      (init_sp*cos(init_co))-x_obs(4,1);];
    P_11 = (init_r^2)*R*(cos(z(1))^2) + (sigma_r^2)*(sin(z(1))^2);
    P_22 = (init_r^2)*R*(sin(z(1))^2) + (sigma_r^2)*(cos(z(1))^2);
    P_12 = ((sigma_r^2) - ((init_r^2)*R))*sin(z(1))*cos(z(1));
    P_21 = P_12;
    P_33 = ((init_sp^2)*(sigma_c^2)*(cos(init_co)^2)) + (sigma_sp^2)*(sin(init_co)^2);
    P_44 = ((init_sp^2)*(sigma_c^2)*(sin(init_co)^2)) + (sigma_sp^2)*(cos(init_co)^2);
    P_34 = ((sigma_sp^2) - (init_sp^2)*(sigma_c^2))*sin(init_co)*cos(init_co);
    P_43 = P_34;
    
    
    P_xx = [P_11,P_12,0,0;...
            P_21,P_22,0,0;...
            0,0,P_33,P_34;...
            0,0,P_43,P_44;];
     
   error_pos(mm,1) = ((x_tar(1,1) - (x_rel_obs(1,1) +  x_obs(1,1)))^2) + ((x_tar(2,1) - (x_rel_obs(2,1) +  x_obs(2,1)))^2);
   
   error_vel(mm,1) = ((x_tar(3,1) - (x_rel_obs(3,1) +  x_obs(3,1)))^2) + ((x_tar(4,1) - (x_rel_obs(4,1) +  x_obs(4,1)))^2);
   
   %**************************************************************************
    % SRF algorithm starts here

    for t = 2:kmax

    %Prediction
        
        x_rel_obs(:,t) = (F*x_rel_obs(:,t-1))  - [x_obs(1,t) - x_obs(1,t-1) - (T*x_obs(3,t-1));...
                                                x_obs(2,t) - x_obs(2,t-1) - (T*x_obs(4,t-1));...
                                                x_obs(3,t) - x_obs(3,t-1);...
                                                x_obs(4,t) - x_obs(4,t-1)];
                                     
       P_xx = (F*P_xx*F') + Q;
       
       Q_m = R*((x_rel_obs(1,t)^2) + (x_rel_obs(2,t)^2) + P_xx(1,1)+ P_xx(2,2)).*eye(2);
                                         
                                   
       V = (H*P_xx*(H')) + Q_m;
       
    % Correction
    
        K = P_xx * (H') * inv(V);
        b = [sin(z(t)); cos(z(t))];
        u = (((b')*inv(V)*b)^(-0.5)) * (b') * inv(V) * (H*x_rel_obs(:,t));
        
    
        rho = (u+ (sqrt(2*pi)*((u^2) + 1)* exp(0.5*(u^2))*normcdf(u)))/...
                    (1+ (sqrt(2*pi)*u*exp(0.5*(u^2))* normcdf(u)));
        
        
        gamma = (((b')*inv(V)*b)^(-0.5)) * rho;
        delta = ((b' * inv(V) * b)^(-1)) * (2 + (u*rho) - (rho^2));
        
        x_rel_obs(:,t) = ((eye(4) - K*H)*x_rel_obs(:,t)) + (gamma*K*b); 
        
        P_xx = ((eye(4) - K*H)*P_xx) + (delta*K*b*(b')*(K'));
        
       % error_xpos(mm,t) = x_tar(1,t) - (x_rel_obs(1,t) +  x_obs(1,t));
       % error_ypos(mm,t) = x_tar(2,t) - (x_rel_obs(2,t) +  x_obs(2,t));
       % error_xdot(mm,t)  = x_tar(3,t) - (x_rel_obs(3,t) +  x_obs(3,t));
       % error_ydot(mm,t)  = x_tar(4,t) - (x_rel_obs(4,t) +  x_obs(4,t));
       
       error_pos(mm,t) = ((x_tar(1,t) - (x_rel_obs(1,t) +  x_obs(1,t)))^2) + ((x_tar(2,t) - (x_rel_obs(2,t) +  x_obs(2,t)))^2);
       error_vel(mm,t) = ((x_tar(3,t) - (x_rel_obs(3,t) +  x_obs(3,t)))^2) + ((x_tar(4,t) - (x_rel_obs(4,t) +  x_obs(4,t)))^2);
    end
    
    if sqrt(error_pos(mm,kmax)) > 1
        fprintf('filter diverged, erros > threshold \n')
        track_loss = track_loss +1;
        error_pos(mm,:) = zeros(1,kmax);
        error_vel(mm,:) = zeros(1,kmax);
    else 
        fprintf('passed \n')
    end
end

x_fil_org = x_rel_obs + x_obs;

error_pos = error_pos(any(error_pos,2),:);
error_vel = error_vel(any(error_vel,2),:);
RMSE_pos = ((1/(MC-track_loss)).*sum(error_pos,1)).^0.5;
RMSE_vel = ((1/(MC-track_loss)).*sum(error_vel,1)).^0.5;

%[theta_obs, rho_obs] = cart2pol(x_obs(1,:),x_obs(2,:));
%[theta_tar, rho_tar] = cart2pol(x_tar(1,:),x_tar(2,:));
%[theta_est,rho_est] = cart2pol(x_fil_org(1,:),x_fil_org(2,:));

fprintf('No of times track diverged is %d of %d MC runs\n', track_loss,MC);
fprintf('Track loss for %d Monte carlo runs:  %.4f %% \n', MC, track_loss*100/MC);

file_out = 'RMSE_results.xlsx';

xlswrite(file_out,{'RMSE Position in kilometers',},sheet,'A1');
xlswrite(file_out,RMSE_pos',sheet,'A2');

xlswrite(file_out,{'RMSE Velocity in kmps'},sheet,'B1');
xlswrite(file_out,RMSE_vel',sheet,'B2');


xlswrite(file_out,{'No. of MC runs'},sheet,'C1');
xlswrite(file_out,MC,sheet,'C2');


xlswrite(file_out,{'Track loss'},sheet, 'D1');
xlswrite(file_out,track_loss,sheet, 'D2');

xlswrite(file_out,{'Track loss %'},sheet, 'E1');
xlswrite(file_out,(track_loss*100/MC),sheet, 'E2');


figure
plot(x_obs(1,:),x_obs(2,:),'b'); % Observer Plot X vs Y
hold on
plot(x_tar(1,:),x_tar(2,:),'r'); % Target Plot X vs Y
hold on
plot(x_fil_org(1,:),x_fil_org(2,:),'g','linewidth',1.5) % Est Target Plot X vs Y
hold on
plot(x_obs(1,1),x_obs(2,1),'r*');
hold on
plot(x_fil_org(1,1),x_fil_org(2,1),'r*');
hold on
plot(x_tar(1,1),x_tar(2,1),'r*');
xlabel('Distance in kms','fontsize',14);
ylabel('Distance in kms','fontsize',14);
%legend('Observer','Target','estimated Target','initial observer position','initial estimated position');
title('Observer and Target ','fontsize',14);

%{
figure
plot(RMSE_pos);
ylabel('position RMSE in kms');
xlabel('time in minutes');
%axis([16 30 0 1.5])
title('position RMSE');

figure
plot(RMSE_vel);
ylabel('velocity RMSE in km/sec ');
xlabel('time in minutes');
%xlim([21 30]);
title('velocity RMSE');

figure
polarplot(theta_obs,rho_obs,'r','linewidth',1.5);
hold on
polarplot(theta_tar,rho_tar,'b','linewidth',1.5);
hold on
polarplot(theta_est,rho_est,'g','linewidth',1.5);
%}