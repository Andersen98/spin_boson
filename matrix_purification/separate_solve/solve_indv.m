function p_t = solve_indv(params,k,c_k,t0,tf)
%SOLVE_INDV Outputs the population difference for a value of k (the
%ancillery oscillator) from t0 to tf

bath_elements = params.bath_elements;
num_bath = prod(bath_elements(3,:));

initial_state = zeros(2*num_bath,1);
initial_state(k,1) = c_k;



 H_sb_w = @(psi) H_sb(bath_elements,psi);
 H_b_w = @(psi) H_b(bath_elements,psi);
 H_s_w = @(psi) H_s(bath_elements, psi, params.epsilon);
 FULL_H = @(t,y) -1i.*(H_sb_w(y) + H_s_w(y) + H_b_w(y));
 
 [t,y] = ode15s(@(t,y) FULL_H(t,y),[t0,tf],initial_state);
 p_t = zeros(length(t),2);
 p_t(:,1) = t;
 for ii = 1:length(t)
    rho_spin = TrX(y(ii,:)',2,[2,num_bath]);
    p_t(ii,2) = trace(rho_spin*[1,0;0,-1]);
end
end

