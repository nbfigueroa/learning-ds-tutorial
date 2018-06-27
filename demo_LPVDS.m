%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Demo Script to learn LPV-DS with different GMM fitting options    %
%      and different optimization variants on Drawn or LASA Datasets     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 1:  Draw with GUI (robot sim enabled) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc

fig1 = figure('Color',[1 1 1]);
% Axis limits
limits = [-6 0.5 -0.5 2];
axis(limits)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.55, 0.2646 0.4358]);
grid on

% Draw Reference Trajectories
[data, hp] = draw_mouse_data_on_DS(fig1, limits);
Data = []; x0_all = []; x0_end = []; Data_sh = [];
for l=1:length(data)    
    data_ = data{l};
    x0_end = [x0_end data_(1:2,end)];
    Data = [Data data_];
    x0_all = [x0_all data_(1:2,1)];
    
    % Shift data to origin for (O2)
    data_(1:2,:) = data_(1:2,:) - repmat(data_(1:2,end), [1 length(data_)]);
    data_(3:4,end) = zeros(2,1);

    Data_sh = [Data_sh data_];
end

% Position/Velocity Trajectories
Xi_ref     = Data(1:2,:);
Xi_dot_ref = Data(3:end,:);

% Global Attractor of DS
att_g = mean(x0_end,2);
scatter(att_g(1),att_g(2),100,[0 0 0],'d'); hold on;

% Position/Velocity Trajectories
delete(hp)
if exist('h_att_g','var');  delete(h_att_g); end
[h_data, h_att, h_vel] = plot_reference_trajectories(Data, att_g, [], 10);
grid on;
box on;
title('Demonstrated Trajectories','Interpreter','LaTex','FontSize',20);
xlabel('$x_1$','Interpreter','LaTex','FontSize',20);
ylabel('$x_2$','Interpreter','LaTex','FontSize',20);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 2: Choose from LASA DATASET %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
clear all; close all; clc
[demos, limits] = load_LASA_dataset();
att_g = [0 0]';
sample = 2;
Data = []; x0_all = [];
for l=1:3   
    % Check where demos end and shift    
    data_ = [demos{l}.pos(:,1:sample:end); demos{l}.vel(:,1:sample:end);];    
    Data = [Data data_];
    x0_all = [x0_all data_(1:2,20)];
    clear data_
end

% Position/Velocity Trajectories
Xi_ref     = Data(1:2,:);
Xi_dot_ref = Data(3:end,:);
figure('Color',[1 1 1])
[h_data, h_att, h_vel] = plot_reference_trajectories(Data, att_g, [], 10);
grid on;
box on;
title('Demonstrated Trajectories','Interpreter','LaTex','FontSize',20);
xlabel('$x_1$','Interpreter','LaTex','FontSize',20);
ylabel('$x_2$','Interpreter','LaTex','FontSize',20);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Step 1: Fit GMM to Trajectory Data        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GMM Estimation Algorithm %%%%
% 1: GMM-EM Model Selection via BIC
% 2: GMM via Competitive-EM
% 3: CRP-GMM (Collapsed Gibbs Sampler)

est_options = [];
est_options.type        = 1;   % GMM Estimation Alorithm Type    
est_options.maxK        = 15;  % Maximum Gaussians for Type 1/2
est_options.do_plots    = 1;   % Plot Estimation Statistics
est_options.adjusts_C   = 1;   % Adjust Sigmas
est_options.fixed_K     = [];  % Fix K and estimate with EM
est_options.exp_scaling = 1;   % Scaling for the similarity to improve locality

% Discover Local Models
sample = 1;
[Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref(:,1:sample:end), Xi_dot_ref(:,1:sample:end), est_options);

%% Extract Cluster Labels
show_pdf = 0;
est_K      = length(Priors0); 
Priors = Priors0; Mu = Mu0; Sigma = Sigma0;
[~, est_labels] =  my_gmm_cluster(Xi_ref, Priors, Mu, Sigma, 'hard', []);

%%% Visualize GMM pdf from learnt parameters (for 1D Datasets)
clear ds_gmm; ds_gmm.Mu = Mu; ds_gmm.Sigma = Sigma; 
ds_gmm.Priors = Priors; 

% Adjust Covariance Matrices
if est_options.adjusts_C  == 1
    tot_scale_fact = 1; rel_scale_fact = 0.25;
    Sigma = adjust_Covariances(Sigma0, tot_scale_fact, rel_scale_fact);
    ds_gmm.Sigma = Sigma;
end    

% Visualize Cluster Parameters on Manifold Data
plotGMMParameters( Xi_ref, est_labels, Mu, Sigma);
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
switch est_options.type
    case 0
        title('Physically-Consistent Non-Parametric Mixture Model','Interpreter','LaTex', 'FontSize',15);
    case 1
        title('Best fit GMM with EM-based BIC Model Selection','Interpreter','LaTex', 'FontSize',15);
    case 3
        title('Bayesian Non-Parametric Mixture Model','Interpreter','LaTex', 'FontSize',15);
end

if show_pdf 
ml_plot_gmm_pdf(Xi_ref, Priors, Mu, Sigma, limits)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 2: ESTIMATE SYSTEM DYNAMICS MATRICES  %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%
% Type of constraints/optimization 
constr_type = 2;      % 0:'convex':     A' + A < 0
                      % 1:'non-convex': A'P + PA < 0
                      % 2:'non-convex': A'P + PA < -Q given P                                  
init_cvx    = 1;      % 0/1: initialize non-cvx problem with cvx                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if constr_type == 0 || constr_type == 1
    P_opt = [];
else
    [Vxf] = learn_wsaqf(Data, att_g);
    P_opt = Vxf.P(:,:,1);
end

%%%%%%%%  LPV system sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%            
if constr_type == 1
    [A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data_sh, [0 0]', constr_type, ds_gmm, P_opt, init_cvx);
    ds_lpv = @(x) lpv_ds(x-att_g, ds_gmm, A_g, b_g);
else
    [A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data, att_g, constr_type, ds_gmm, P_opt, init_cvx);
    ds_lpv = @(x) lpv_ds(x, ds_gmm, A_g, b_g);
end


%% %%%%%%%%%%%%    Plot Resulting DS  %%%%%%%%%%%%%%%%%%%
fig1 = figure('Color',[1 1 1]);
[hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 0 0],'filled'); hold on
h_ds = plot_ds_model(fig1, ds_lpv, [0 0]', limits,'medium'); hold on;
scatter(att_g(1),att_g(2),150,[0 0 0],'d', 'MarkerFaceColor',[1 1 1]); hold on;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on
title('LPV-DS $\mathbf{f}(x)$', 'Interpreter','LaTex','FontSize',20)
xlabel('$x_1$','Interpreter','LaTex','FontSize',20);
ylabel('$x_2$','Interpreter','LaTex','FontSize',20);

% Simulate trajectories and plot them on top
plot_repr = 1;
if plot_repr
    opt_sim = [];
    opt_sim.dt = 0.01;
    opt_sim.i_max = 3000;
    opt_sim.tol = 0.1;
    opt_sim.plot = 0;
    [x_repr, ~]=Simulation(x0_all ,[],ds_lpv, opt_sim);
        if exist('hr','var');     delete(hr);    end
    [hr] = scatter(x_repr(1,:),x_repr(2,:),10,[0 0 0],'filled'); hold on
    
    % Compute DTWD between train trajectories and reproductions
    nb_traj       = size(x_repr,3);
    ref_traj_leng = size(Xi_ref,2)/nb_traj;
    dtwd = zeros(1,nb_traj);
    for n=1:nb_traj
        start_id = 1+(n-1)*ref_traj_leng;
        end_id   = n*ref_traj_leng;
        dtwd(1,n) = dtw(x_repr(:,:,n)',Xi_ref(:,start_id:end_id)',20);
    end
    fprintf('LPV-DS got DTWD of reproduced trajectories: %2.4f +/- %2.4f \n', mean(dtwd),std(dtwd));
    
end

% Compute RMSE on training data
rmse = mean(rmse_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got prediction RMSE on training set: %d \n', constr_type+1, rmse);

% Compute e_dot on training data
edot = mean(edot_error(ds_lpv, Xi_ref, Xi_dot_ref));
fprintf('LPV-DS with (O%d), got e_dot on training set: %d \n', constr_type+1, edot);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Plot Choosen Lyapunov Function and derivative  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of plot
contour = 1; % 0: surf, 1: contour
clear lyap_fun_comb lyap_der 

switch constr_type
    case 0 
        P = eye(2);
        title_string = {'$V(x) = (x-x^*)^T(x-x^*)$'};
    case 1
        P = P_est;
        title_string = {'$V(x) = (x-x^*)^TP(x-x^*)$'};
    case 2
        P = P_opt;
        title_string = {'$V(x) = (x-x^*)^TP(x-x^*)$'};
end

% Lyapunov function
lyap_fun = @(x)lyapunov_function_PQLF(x, att_g, P);

% Derivative of Lyapunov function (gradV*f(x))
lyap_der = @(x)lyapunov_derivative_PQLF(x, att_g, P, ds_lpv);
title_string_der = {'Lyapunov Function Derivative $\dot{V}(x)$'};

if exist('h_lyap','var');     delete(h_lyap);     end
if exist('h_lyap_der','var'); delete(h_lyap_der); end
h_lyap     = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
[h_data, h_att, h_vel] = plot_reference_trajectories(Data, att_g, [], 10);
xlabel('$x_1$','Interpreter','LaTex','FontSize',20);
ylabel('$x_2$','Interpreter','LaTex','FontSize',20);
h_lyap_der = plot_lyap_fct(lyap_der, contour, limits_,  title_string_der, 1);
[hd] = scatter(Xi_ref(1,:),Xi_ref(2,:),10,[1 1 0],'filled'); hold on
xlabel('$x_1$','Interpreter','LaTex','FontSize',20);
ylabel('$x_2$','Interpreter','LaTex','FontSize',20);

%% Compare Velocities from Demonstration vs DS
% Simulated velocities of DS converging to target from starting point
xd_dot = []; xd = [];
% Simulate velocities from same reference trajectory
for i=1:length(Data)
    xd_dot_ = ds_lpv(Data(1:2,i));    
    % Record Trajectories
    xd_dot = [xd_dot xd_dot_];        
end

% Plot Demonstrated Velocities vs Generated Velocities
if exist('h_vel','var');     delete(h_vel);    end
h_vel = figure('Color',[1 1 1]);
plot(Data(3,:)', '.-','Color',[0 0 1], 'LineWidth',2); hold on;
plot(Data(4,:)', '.-','Color',[1 0 0], 'LineWidth',2); hold on;
plot(xd_dot(1,:)','--','Color',[0 0 1], 'LineWidth', 1); hold on;
plot(xd_dot(2,:)','--','Color',[1 0 0], 'LineWidth', 1); hold on;
grid on;
legend({'$\dot{x}^{ref}_{1}$','$\dot{x}^{ref}_{2}$','$\dot{x}^{d}_{1}$','$\dot{x}^{d}_{2}$'}, 'Interpreter', 'LaTex', 'FontSize', 15)
