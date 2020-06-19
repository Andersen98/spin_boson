clear variables;
%energy is in units of Delta so we are left with the
%physics params - (k_b=h_bar=1, energy in units of delta)
params.epsilon = 0;
params.lambda = 2.5;
params.T = 2; 

%simulation parameters
    %osc_count is the number of oscillators to be used in the bath
    %n_max is the maximal mode for each oscillator i,e, can range from
    %0,1,2,3...n
params.osc_count = 6;
params.mode_counts = 3.*ones(1,params.osc_count);

%Descretized Debey spectral density - need us to choose a
    %charactaristic frequency for the bath wc
    %core frequency range: min to max [w0,wf]
    %coupling strength (renormalizatin energy) lambda
params.wc= 4;
params.wf = 2*params.wc;
params.w0 = params.wf./(2.*params.osc_count);
params.lambda = 2.5;

%construct the bath terms (frequency and coupling
bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.mode_counts,params.lambda);

%get a truncated initial bath state, where the bath is in "thermal
%equilibrium"


%make an initial state that would give a initial density matrix \rho(0) in
%section IV of Berkelbach 2011
disp("making initial state")
tic()
initial_state = psi_naught(bath_elements,1/params.T,[1;0]);
toc()

disp("create hamiltonian")
tic()
H_couple = H_sb(bath_elements);
H_spin = H_s(bath_elements,params.epsilon);
t_span = [0,12];
toc()

disp("solve ode")
tic()
H = H_spin;
[t,y] = ode45(@(t,y) -1i.*(H*y),t_span,initial_state);
toc()

disp("take partial trace")
tic()
%helper variables
osc_count = length(bath_elements(1,:));
n_bath = prod(params.mode_counts);
n_total = 2*n_bath;
pt = zeros(length(t),1);
sys = [2,params.mode_counts];
for ii = 1:length(t)
    rho_spin = TrX(y(ii,:)',2:(osc_count+1),sys);
    %rho_spin = [sum(y(ii,1:n_bath).^2),0;0,sum(y(ii,n_bath+1:n_total).^2)];
    pt(ii) = trace(rho_spin*[1,0;0,-1]);
end
toc()
%plot the solution
plot(t,pt(:,1),'-o');
title('Poulation difference Tr[\rho_s \sigma_z] solved with ODE 45')
xlabel('Time t');
ylabel('Population Difference');
