clear variables;
%energy is in units of Delta so we are left with the
%physics params - (k_b=h_bar=1, energy in units of delta)
params.epsilon = 0;
params.lambda = 2.5;
params.T = .2; 

%simulation parameters
    %osc_count is the number of oscillators to be used in the bath
    %n_max is the maximal mode for each oscillator i,e, can range from
    %0,1,2,3...n
params.n_max = 1;
params.osc_count = 10;

%Descretized Debey spectral density - need us to choose a
    %charactaristic frequency for the bath wc
    %core frequency range: min to max [w0,wf]
    %coupling strength (renormalizatin energy) lambda
params.wc= 4;
params.wf = 2*params.wc;
params.w0 = params.wf./(2.*params.osc_count);
params.lambda = 2.5;

%construct the bath terms (frequency and coupling
bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.lambda);

%spin coupling Hamiltonian
% H_debye(params.n_max,bath_elements);
%bath Hamiltonian
% H + H_b(params.n_max,bath_elements);
%spin Hamiltonian
% H + H_s(params.n_max,bath_elements,params.epsilon);
H =  H_s(params.n_max,bath_elements,params.epsilon) + H_b(params.n_max,bath_elements) + H_debye(params.n_max,bath_elements);

%make an initial state that would give a initial density matrix \rho(0) in
%section IV of Berkelbach 2011
initial_state = psi_naught(params.n_max,bath_elements,1/params.T);

%helper variables
mode_count = params.n_max +1;
osc_count = length(bath_elements(1,:));
n_bath = mode_count^osc_count;
n_total = 2*n_bath;

%now evolve state
n_steps =150;
t = linspace(0,12,n_steps);
p_diff = zeros(1,n_steps);
tic()
parfor ii = 1:n_steps
    psi_ii = expm(-1i.*t(ii).*H)*initial_state;
    %psi_up_ii = psi_ii(1:n_bath);
   % psi_down_ii = psi_ii((n_bath+1):n_total);
    %p_diff(ii) =  norm(psi_ii(1:n_bath))^2 - norm(psi_ii((n_bath+1):n_total))^2 ;
    p_diff(ii) = trace(TrX(psi_ii,2,[2,n_bath])*[1,0;0,-1]);
    psi_ii = [];
end
toc()
plot(t,p_diff)
    