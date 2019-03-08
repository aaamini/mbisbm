%% Set-up
addpath('../common')
clear all
K = 10;
N = 2*[100,400]; %N should be devisable by K

sig2 = 5*[.1,.1];

d = [2 2];
nu = 10;
mu = zeros(1,sum(d));
TAG = '_v2_A';


mylogspace = @(x,y,n) logspace(log10(x),log10(y),n);

lam_len = 12;
lam_vec = linspace(2,40,lam_len);

beta = 1/7;
nMC = 5;

[nmi_cell, miss_cell] = deal( cell(lam_len,nMC) );
deg_vec = zeros(lam_len,1);   

random_init = @(N,K) mnrnd(1,ones(1,K)/K,N);
rho = 0.1;

methods = {'~rnd','biSC','SCP', ...
    'mbiSBM (~rnd)','mbiSBM (biSC)','mbiSBM (SCP)', ...
    'mbiSBM (~rnd)','mbiSBM (biSC)','mbiSBM (SCP)'
    };


%%
for r = 1:lam_len
    fprintf('--- r = %3d ---\n', r)
    lambda = lam_vec(r);
    
    temp_deg = 0;
    parfor tmc = 1:nMC
         while 1 % ~SUCCESS
            try
                if mod(tmc,10) == 0, fprintf('\n'), end
                fprintf('.')

                Sig = nu*eye(sum(d));

                %[ A, label_1, label_2, X_1, X_2, v ] = bip_sbm_mix2( N, K, p, q, pr, d, mu ,Sig, sig2 );
                %[A,label_1,label_2,theta_1,theta_2,expected_A] = bip_dcsbm(N, K, p, q, pr, delta);
                %[A, label, theta, expected_A] = bip_dcsbm2(N, K, p, q, pr, 'theta_low', .1, 'theta_high', .5, 'alpha', .1);
                [A, label, theta, emp_avg_deg, Q, pr, expected_A] = bip_dcsbm_pareto(N, K, beta, lambda, 'alpha', inf, 'poisson', false);
                v = zeros(K,d(1)+d(2));
                for k=1:K
                    v(k,:) = mvnrnd(mu,Sig);
                end
                idx = 1:d(1);
                X_1 = v(label{1},idx) +  sig2(1)*randn(N(1),d(1));
                idx = d(1) + (1:d(2));
                X_2 = v(label{2},idx) +  sig2(2)*randn(N(2),d(2));

                tru_lab = BiClustRes.setget_tru_lab(label);  % label is a 1x2 cell array 
                BiClustRes.setget_comp_miss(false);

                [z1_rnd, z2_rnd] = dirchlet_perturb(label_vec2mat(tru_lab{1}) , label_vec2mat(tru_lab{2}) , 1-rho, .5);

                temp_deg = temp_deg + emp_avg_deg;

                Xs = {X_1;X_2};
                res = {};

                %% ~ rnd
                res{end+1} = BiClustRes({z1_rnd; z2_rnd},'~rnd');

                %%% biSC
                NORMALIZE = 1;
                [l1_bisc, l2_bisc, Z_2] = biSpecClust(A,K,NORMALIZE,1:K);
                res{end+1} = BiClustRes({l1_bisc; l2_bisc},'biSC');
                z1_bisc = label_vec2mat(l1_bisc);
                z2_bisc = label_vec2mat(l2_bisc);

                %%% SCP
                At = [zeros(N(1),N(1)) A;A' zeros(N(2),N(2))];        
                e = initLabel5b(At, K, 'scp', struct('verbose',false));
                l1_scp = e(1:N(1));
                l2_scp = e(N(1) + (1:N(2)));
                res{end+1} = BiClustRes({l1_scp; l2_scp},'SCP');
                z1_scp = label_vec2mat(l1_scp,K);
                z2_scp = label_vec2mat(l2_scp,K);

                %%  mbiSBM (no X)
                opts = {'verbose',1,'alpha',.1,'mu',0.01,'first_tau_update',1, 'init_p_q_pr', true};
                [tau_1, tau_2] = fit_mbiSBM_v2(A, {[];[]}, K, z1_rnd, z2_rnd, 'ignore_theta',1, opts{:});
                res{end+1} = BiClustRes({tau_1; tau_2},'mbiSBM (~rnd), no X');

                [tau_1, tau_2] = fit_mbiSBM_v2(A, {[];[]}, K, z1_bisc, z2_bisc, 'ignore_theta',1, opts{:});
                res{end+1} = BiClustRes({tau_1; tau_2},'mbiSBM (biSC), no X');

                [tau_1, tau_2] = fit_mbiSBM_v2(A, {[];[]}, K, z1_scp, z2_scp, 'ignore_theta',1, opts{:});
                res{end+1} = BiClustRes({tau_1; tau_2},'mbiSBM (SCP), no X');

                %%  mbiSBM (X)
                opts = {'verbose',1,'alpha',.1,'mu',0.01,'first_tau_update',1, 'init_p_q_pr', true};
                [tau_1, tau_2] = fit_mbiSBM_v2(A, Xs, K, z1_rnd, z2_rnd, 'ignore_theta',1, opts{:});
                res{end+1} = BiClustRes({tau_1; tau_2},'mbiSBM (~rnd), X');

                [tau_1, tau_2] = fit_mbiSBM_v2(A, Xs, K, z1_bisc, z2_bisc, 'ignore_theta',1, opts{:});
                res{end+1} = BiClustRes({tau_1; tau_2},'mbiSBM (biSC), X');

                [tau_1, tau_2] = fit_mbiSBM_v2(A, Xs, K, z1_scp, z2_scp, 'ignore_theta',1, opts{:});
                res{end+1} = BiClustRes({tau_1; tau_2},'mbiSBM (SCP), X');


                nMTD = length(res); 
                [nmi_temp, miss_temp] = deal(zeros(nMTD,1));
                for k = 1:nMTD
                    nmi_temp(k) = res{k}.nmi; 
         %           miss_temp(k) = res{k}.miss_rate;
                end

                nmi_cell{r,tmc} = nmi_temp; 
        %        miss_cell{r,tmc} = miss_temp; 
                break
            catch
                % do nothing
                disp('---------------- error occued, skipping ---------------------')
            end
         end
    end
    fprintf('\n')
    
    deg_vec(r) = temp_deg / (nMC);
end

%%
nMTD = length(nmi_cell{1,1});
myreshape = @(X) shiftdim( reshape([X{:}],[nMTD size(X)]), 1);
nmi = myreshape(nmi_cell);
%miss = myreshape(miss_cell);

%save(['NMI_deg_covar_nu_sepNMI_' TAG])
%file_name = ['../results/compare_nmi_bm' TAG];
%save(file_name)

%%
figure(3), clf, 
set(gcf, 'PaperPositionMode','manual')
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 7 5.5])

markers = {'o-','^-','*--','s-','d--','x-.'};
cc = jet(6);
 for k = 1:6
    h(k) = plot(deg_vec,mean(nmi(:,:,k),2),markers{k});
    hold on
end
legend(methods,'location','southeast','fontsize',12);
legend('location','southeast')
legend('boxoff')
ylabel('Matched NMI'), xlabel('Average degree'), axis tight, ylim([0,1]),
set(h,'LineWidth',1.0)

plot_multi_line_text({  sprintf('$d = (%d,%d)$',d(1),d(2)); 
                        sprintf('$K = %d$',K);
                        sprintf('$\\alpha = %3.2f$', beta);
                        sprintf('$\\nu = %3.1f$',0)},'left',.15,'bottom',.7)


%print('-depsc',[file_name '_no_X.eps'])                    

%%
figure(4), clf, 
set(gcf, 'PaperPositionMode','manual')
set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 7 5.5])

markers = {'o-','^-','*--','s-','d--','x-.'};
cc = jet(6);
mthd_idx = [1:3 7:9];
 for k = 1:6
    h(k) = plot(deg_vec,mean(nmi(:,:,mthd_idx(k)),2),markers{k});%, 'color',cc(k,:))
    hold on
end
legend(methods,'location','southeast','fontsize',12);
legend('location','southeast')
legend('boxoff')
ylabel('Matched NMI'), xlabel('Average degree'), axis tight,  ylim([0,1])
set(h,'LineWidth',1.0)

plot_multi_line_text({  sprintf('$d = (%d,%d)$',d(1),d(2)); 
                        sprintf('$K = %d$',K);
                        sprintf('$\\alpha = %3.2f$', beta);
                        sprintf('$\\nu = %3.1f$',nu)},'left',.4,'bottom',.15)

%print('-depsc',[file_name '_has_X.eps'])



