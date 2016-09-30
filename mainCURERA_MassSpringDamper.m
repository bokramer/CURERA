% Main file for the numerical experiments in the paper
% "System Identifcation via CUR-factored Hankel approximation", 2016, by
% Boris Kramer and Alex A. Gorodetsky 
%
%
%
% Copyright (c) MIT, 2016
% Boris Kramer (bokramer@mit.edu)
% ------------------------------------------------------------------------

clear

load('MSD_model','A','B','C')     % m=p=30, N=1000, sparse, no D term
format shortE
savefigs = 0;       % =1 if you want figures to be saved 

%% Initialization
 m = size(B,2);     % Inputs
 p = size(C,1);     % Outputs
 N = max(size(A));  % State space dimension
 T = 2*1e4;         % 
Ns = 500;           % Collecting 2*Ns samples though
Ts = T/(2*Ns);      % Sampling time
 r = 80;            % Reduced-order model dimension

%% % Discrete Time System for sampling of Markov Parameters
sysfull  = ss(full(A),full(B),full(C),0);
sysfulld = c2d(sysfull,Ts,'tustin');        % bilinear transformation
[Ad,Bd,Cd,Dd] = ssdata(sysfulld);           % Now a discrete time system

markov = cell(1,2*Ns);
markovnorm = zeros(1,2*Ns);
markov{1} = Dd;         % First Markov parameter is D term
markovnorm(1) = norm(Dd);

f = Bd;
for jj = 2:2*Ns         % Computing Markov Parameters iteratively 
    g = Cd*f;
    f = Ad*f;           % markov{jj} = Cd*Ad^(jj-1)*Bd;
    markov{jj} = g;     % + 0.02*norm(markov{ii})*randn(p,m); for noisy data
    markovnorm(jj) = norm(g,'fro')^2;
end


%% SVD-based ERA for system ID -- for comparison purposes
H = block_Hankel(markov(2:Ns+1),markov(Ns+1:2*Ns));
Hnorm = norm(H,'fro');   % Frobenius norm of Hankel matrix needed later 

load('svdMSDred','SigZ','U','V')
S = SigZ;
%[U,S,V] = svd(H,'econ');
[Ard,Brd,Crd] = era(U,S,V,r,m,p); 

% Computing terror in the SVD approximation (optimal error)
S = diag(S);
relerrH_SVD = full(sqrt(sum(S(r+1:end).^2))/Hnorm); % ||H-Happrox||_F/||H||F
display(relerrH_SVD)

markoverrorEra = zeros(1,2*Ns);
f=Brd;
for jj = 2:2*Ns          % Collecting Markov Parameters, iterative computation  
    g = Crd*f;
    f = Ard*f;            % markov{i} = Cd*Ad^(i-1)*Bd;
    markoverrorEra(jj) = norm(g-markov{jj},'fro')^2;
end

% Computing error in approximating the Markov sequence
MarkErrERA = sum(markoverrorEra)/sum(markovnorm);
display(MarkErrERA)


%% Low rank matrix approximation through CUR
epsTol = 2e-2;              % Tolerance for maxvol algorithm
rr = r; 
itCUR=1;                    % Since initial column/row selection in maxvol is random,   
tCUR = zeros(itCUR,1);      % can run multiple experiments
Hdiff_full = cell(itCUR,1);
condH = cell(itCUR,1);
Hdiff_iterates = cell(itCUR,1);
HnormIterates = cell(itCUR,1);
relerrH_cur = zeros(itCUR,1);

relTolDelta = 1e-4; % ||(CUR)_{i} - (CUR)_{i+1} ||_F/||(CUR)_{i}||_F < relTolDelta to stop iteration
compMore =0;
for ii =1:itCUR
    t1 =cputime;
    [I,J,u,s,v, Hdiff_full{ii,jj}, condH{ii,jj}, Hdiff_iterates{ii,jj}, HnormIterates{ii,jj}] = crossApprox(H,rr,relTolDelta,epsTol, compMore);
    t2 = cputime;
    tCUR(ii) = t2-t1;
    
    Hnew = ((H(:,J)*v)*diag((1./diag(s)))*((u')*H(I,:)));
    relerrH_cur(ii) = norm(H-Hnew,'fro')/Hnorm;
end

% Computational time assessment
TCUR = sum(tCUR)/itCUR;
display(TCUR)

% Relative error ||(CUR) - H||_F/||H||_F over all iteration
display(sum(relerrH_cur)/itCUR);


%% CUR-based ERA
M1 = H(:,J)*v;
M2 = (u')*H(I,:);
[qO,rO] = qr(M1,0);
[qC,rC] = qr(M2',0);
[u1,s1,v1] = svd(rO*(diag(1./diag(s)))*(rC'));
Uhat = qO*u1;
Vhat = qC*v1;

[Ard1,Brd1,Crd1] = era(Uhat,s1,Vhat,r,m,p); 


%% Computing Markov errors
markoverrorLR = zeros(1,2*Ns);
f=Brd1;
for jj = 2:2*Ns          % Collecting Markov Parameters, iterative computation  
    g = Crd1*f;
    f = Ard1*f;            % markov{i} = Cd*Ad^(i-1)*Bd;  
    markoverrorLR(jj) = norm(g-markov{jj},'fro')^2;
end

% Computing the error in approximating the Markov sequence
MarkerrLR = sum(markoverrorLR)/sum(markovnorm);
display(MarkerrLR);


%% Closer look at the spectrum: Computing spectrum of identified matrix
Lera = eig(Ard);
Lcur = eig(Ard1); 

% Plotting Eigenvalue accuracy
[~,I1] = sort(abs(Lera),'descend');
[~,I2] = sort(abs(Lcur),'descend');
ldiff = abs(abs(Lera(I1)) - (abs(Lcur(I2))))./abs(Lera(I1));

% Eigenvalue Plots
h4 = figure;
theta = linspace(0,2*pi,1000); % unit circle for comparison
plot(real(Lcur),imag(Lcur),'+');  hold on;
plot(sin(theta),cos(theta),'k');
plot(real(Lera),imag(Lera),'ro'); hold off
xlabel('Re'); ylabel('Im')
legend('CUR','SVD','Location','best')
%title('Spectrum of Ard via CUR approximation');
if savefigs == 1
    export_fig -transparent msd_EVlocations.pdf
    saveas(h4, './figures/msd_EVlocations.fig','fig')
    printSpec = '%3.12f %3.12f %3.12f %3.12f\n';
    fid = fopen('eigenvalues.txt', 'w');
    fprintf(fid,'RERA IERA RCUR ICUR\n');
    fprintf(fid,printSpec, [real(Lera) imag(Lera) real(Lcur) imag(Lcur)]);
    fclose(fid);
end

h5 = figure;
loglog(1-abs(Lcur(I2)),ldiff,'+');
xlabel('|1-|\lambda_i^{CUR} |  |');
ylabel('| (\lambda_i^{CUR} - \lambda_i^{ERA})/\lambda_i^{CUR}  |');

if savefigs == 1
    export_fig -transparent msd_EVaccuracy.pdf
    saveas(h5, './figures/msd_Evaccuracy.fig','fig')
    
    printSpec = '%3.12f %3.12f\n';
    fid = fopen('eigenDistUnitCircle.txt', 'w');
    fprintf(fid,'dist ldiff\n');
    fprintf(fid,printSpec, [1-abs(Lcur(I2)) ldiff]);
    fclose(fid);
end


%% Post Processing - Comparison of Transfer functions
sysredera = d2c(ss(Ard,Brd,Crd,Dd,Ts),'tustin'); % mapping back to cont time
[Ar,Br,Cr,~] = ssdata(sysredera);

sysredcur = d2c(ss(Ard1,Brd1,Crd1,Dd,Ts),'tustin');  % and adding the D term back into the model
[Ar1,Br1,Cr1,~] = ssdata(sysredcur);

w = logspace(-4,1,1000);                % frequency range for Bode plot
H = @(s) C*((s*speye(size(A))-A)\B);    % FOM transfer function
Hcur = @(s) Cr1*( (s*speye(size(Ar1)) -Ar1 )\Br1) ; % ROM transfer functions
Hera = @(s) Cr*( (s*speye(size(Ar)) -Ar )\Br);

% Allocation
tfeval = cell(size(w));
norm_tfeval = zeros(size(w));
tfeval_era = cell(size(w));
tfeval_LR = cell(size(w));
norm_tfeval_era = zeros(size(w));
norm_tfeval_CUR = zeros(size(w));

% COmputing magnitude plot of transfer function
for jj = 1:length(w)
    tfeval_era{jj} = full(Hera(1i*w(jj)));
    tfeval_LR{jj} = full(Hcur(1i*w(jj)));
    tfeval{jj} = full(H(1i*w(jj)) );
    norm_tfeval(jj) = norm(tfeval{jj},2);
    norm_tfeval_era(jj) = norm(tfeval_era{jj},2);
    norm_tfeval_CUR(jj) = norm(tfeval_LR{jj},2);
end

h7 = figure;
loglog(w,norm_tfeval,w,norm_tfeval_era,'k--',w,norm_tfeval_CUR,'g.-');
set(gca,'DefaultTextFontSize',18)
xlabel('Frequency w')
legend('||H(iw)||','||H_{era}(iw)||','||H_{CUR}(iw)||','Location','SouthWest')


if savefigs == 1
    saveas(h7,'./figures/msd_TF_all.fig','fig')
    saveas(h7,'./figures/msd_TF_all.eps','epsc2')
    export_fig -transparent msd_msd_TF_all.pdf
    
    printSpec = '%3.12f %3.12f %3.12f %3.12f \n';
    fid = fopen('transferFunction.txt', 'w');
    fprintf(fid,'freq norm_tf_eval norm_tfeval_era norm_tfeval_CUR \n');
    fprintf(fid,printSpec, [w norm_tfeval norm_tfeval_era norm_tfeval_CUR]);
    fclose(fid);
end


%% Compute the relative H2 error
hinferr_era = zeros(size(w));
hinferr_cur = zeros(size(w));
for i = 1: length(w)
    hinferr_era(i) = norm(tfeval_era{i} - tfeval{i},2)/norm(tfeval{i},2);
    hinferr_cur(i) = norm(tfeval_LR{i} - tfeval{i},2)/norm(tfeval{i},2);
end

h8 = figure; set(gca,'fontsize',14)
loglog(w,hinferr_era,w,hinferr_cur,'Linewidth',1.1)
legend('ERA vs Full','CUR vs Full')
xlabel('Frequency w ')
%title('||H(iw) - Hr(iw)||')

if savefigs == 1
    saveas(h8,'./figures/msd_TF_err.fig','fig')
    saveas(h8,'./figures/msd_TF_err.eps','epsc2')
    export_fig -transparent msd_msd_TF_err.pdf
    
    printSpec = '%3.12f %3.12f %3.12f \n';
    fid = fopen('transferFunctionError.txt', 'w');
    fprintf(fid,'freq hinferr_era hinferr_cur \n');
    fprintf(fid,printSpec, [w hinferr_era hinferr_cur]);
    fclose(fid);
end


%%  Continuous time systems
sysredfull = ss(Ad,Bd,Cd,Dd);
[Ad,Bd,Cd,Dd] = ssdata(sysfull);

sysredera = d2c(ss(Ard,Brd,Crd,Dd,Ts),'tustin');
[Ar,Br,Cr,Dr] = ssdata(sysredera);

sysredcur = d2c(ss(Ard1,Brd1,Crd1,Dd,Ts),'tustin');    % adding the D term back into the model
[Ar1,Br1,Cr1,Dr1] = ssdata(sysredcur);

tspan = [0 50];
inp = @ (t) exp(-0.1*t).*sin(5*t); % In paper: exp(-0.05*t).*sin(5*t);
p1 = 1:30; % Number of inputs for the simulation. 


%% Simulate the reduced-order identified systems, compare to FOM
[T_full,Y_full]  = ode45(@(t,y) A*y + B(:,p1)*((inp(t)*ones(length(p1),1))), tspan, zeros(N,1));
Y_out_full = C*Y_full';

[~,Y_era]  = ode45(@(t,y) Ar*y + Br(:,p1)*(inp(t)*ones(length(p1),1)), T_full, zeros(r,1));
Y_out_era = Cr*Y_era';

[~,Y_cur]  = ode45(@(t,y) Ar1*y + Br1(:,p1)*(inp(t)*ones(length(p1),1)), T_full, zeros(r,1));
Y_out_cur = Cr1*Y_cur';

for jj = 1 : 5 : p
    h = figure;  set(gca,'fontsize',18)
    plot( Y_out_full(jj,:),'k','Linewidth',1.05); hold on
    plot(Y_out_era(jj,:),'g--','Linewidth',1.05); hold on;
    plot(Y_out_cur(jj,:),'r-.','Linewidth',1.05)
    legend('Full','ERA','CUR')
    hold off
    xlabel('Time (s)')
    ylabel(strcat('Output No.',num2str(jj) ) );
    if savefigs == 1
        saveas(h,strcat('./figures/msd_CTr80_out',num2str(jj),'noise','.fig'),'fig')
        saveas(h,strcat('./figures/msd_CTr80_out',num2str(jj),'noise','.eps'),'epsc2')
        str = strcat('./figures/msd_CTr80_out',num2str(jj),'.pdf');
        export_fig -transparent string
    end
end

printSpec = '%3.12f %3.12f %3.12f %3.12f \n';
fid = fopen('output6sim.txt', 'w');
fprintf(fid,'T out6fom out6era out6cur \n');
fprintf(fid, printSpec, [T_full Y_out_full(6,:)' Y_out_era(6,:)' Y_out_cur(6,:)']);
fclose(fid);

printSpec = '%3.12f %3.12f %3.12f %3.12f \n';
fid = fopen('output11sim.txt', 'w');
fprintf(fid,'T out6fom out6era out6cur \n');
fprintf(fid, printSpec, [T_full Y_out_full(11,:)' Y_out_era(11,:)' Y_out_cur(11,:)']);
fclose(fid);



