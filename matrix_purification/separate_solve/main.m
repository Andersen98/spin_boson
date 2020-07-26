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
params.osc_count =7;
params.mode_count = 2;

%Descretized Debey spectral density - need us to choose a
    %charactaristic frequency for the bath wc
    %core frequency rang    e: min to max [w0,wf]
    %coupling strength (renormalizatin energy) lambda
params.wc= .25;
params.wf = 2*params.wc;
params.w0 = params.wf./(2.*params.osc_count);
params.lambda = 2.5;

%construct the bath terms (frequency and coupling
params.bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.lambda,params.mode_count);
num_bath = prod(params.bath_elements(3,:));
initial_state_terms = psi_naught_terms(params.bath_elements,1/params.T);

t0 = 0;
tf = 12;
n_pts = 150;
tq= linspace(t0,tf,n_pts);
pq = zeros(num_bath,n_pts);
for k = 1:num_bath
    k
    pk = solve_indv(params,k,initial_state_terms(k),t0,tf);
    pq(k,:) = interp1(pk(:,1),pk(:,2),tq,'spline');
    plot(tq,sum(pq,1),'o')
    drawnow
end
