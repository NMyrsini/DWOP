function [obj_u,obj_v,pred_mean,pred_var, nu_next,sigma_stack] = DWOP_update(event,obj_u,obj_v,part,nu,one_,sigma_stack,n,initBool)
%This function takes one event and updates the pair of approximate distributions, applies accept-reject sampling
%and Newton Raphson updates.





%one_ = 1 and it directly passed as an argument in sliding_windows function
dim = length(obj_u.mu);
%initialize the main diagonal of cov. matrices
V_u = one_*eye(length(obj_u.mu))/dim;
V_v = one_*eye(length(obj_v.mu))/dim;
%nu2 = nu^2;

target_sigma = 1.76;
%left boundary of the ordered probit model
left = -10.8;

% prior mean vector (user) at the current time step
prior_mu1 = obj_u.mu;
% prior covariance matrix (user) at the current time step
prior_Sig1 = obj_u.Sig + (event.time - obj_u.time)*V_u;
% prior mean vector (hotel) at the current time step
prior_mu2 = obj_v.mu;
% prior covariance matrix (hotel) at the current time step
prior_Sig2 = obj_v.Sig + (event.time - obj_v.time)*V_v;
%conduct a rating prediction before an observation occurs.
pred_mean = prior_mu1'*prior_mu2;

%this is not used in the current work. It remains for experimantation purposes
pred_var = prior_mu1'*prior_Sig2*prior_mu1 + prior_mu2'*prior_Sig1*prior_mu2 + sum(sum(prior_Sig1.*prior_Sig2));


%posterior moments initialization
post_mu1 = prior_mu1;%user post. mean vector
post_Sig1 = 0*prior_Sig1;%user post. covariance matrix
post_mu2 = prior_mu2;%hotel post. mean vector
post_Sig2 = 0*prior_Sig2;%hotel post. covariance matrix

%Temporary Buffers
post_mu3 = post_mu1;
post_Sig3 = 0*prior_Sig1;
post_mu4 = prior_mu2;
post_Sig4 = 0*prior_Sig2;

% parameters used to calculate Eqs. 49-52 in the paper
invSig1 = pinv(prior_Sig1);
invSig2 = pinv(prior_Sig2);
v1 = invSig1*prior_mu1;
v2 = invSig2*prior_mu2;

% UPDATE distributions
sigma_0 = 1.1;
x0=sigma_0;
flag = 0;
flag2 =0;
%When there is not prior information about a hotel use the segments  
init_mo = (part(1)+part(2))/2;

%initialize values when a hotel appears for the first time
if initBool == 1
    pred_mean = init_mo;
    prev_sig =  target_sigma;
    flag = 1;
else
%Otherwise find the previous std. deviation related to the current hotel under investigation  
    for j = 2 : n-1 
        if sigma_stack(n).HotelName == sigma_stack(n-j).HotelName
            prev_sig = sigma_stack(n-j).sigma; 
            break 
        else
            prev_sig = target_sigma;
        end  
    end

end

nu = prev_sig;
nu2 = prev_sig^2;
exodos =0;
for x=1:5  
    %because in first appearance of hotel j prior_mu1'*prior_mu2 is very
    %high. Initialize once to solve it.
    m = post_mu1'*post_mu2;
    
    if flag==1 || m<left|| m>10
        m=init_mo;
    else
        m = post_mu1'*post_mu2;
    end
	%compute the bounds based on the left  (part(1)) and the right (part(2)) ones which define the segment 
	%w.r.t. a hotel rating (Eq. 29 in the paper)
    alpha = ((part(1) - m)/nu);
    beta = ((part(2) - m)/nu);
	%pdf od std. normal distribution
    phi_a = normpdf(alpha,0,1);
    phi_b = normpdf(beta,0,1);
	%cdf od std. normal distribution
    Phi_a = normcdf(alpha,0,1);
    Phi_b = normcdf(beta,0,1);
 
    
    Z = Phi_b - Phi_a;
    
    %accept - reject sampling scheme to compute the optimal expected value 
    EXX=trandn(alpha,beta); 
    Zi=m+nu*EXX;%Eq. 35 in the paper
    Ey=Zi; 
    
    if isnan(Ey) || isinf(Ey) || Ey < part(1) || Ey > part(2)
        Ey = .1*sign(part(1)+part(2));
    end

    

 %----- optimal sigma through Newton_Raphson Updates      
   x0 = target_sigma;%initial guess of std. deviation
 %----------------------------------------------------  
     if exodos == 0
         %10 iterations ensure convergence of NR update. 
         for i=1:10
             cdf1 = ((phi_a - phi_b)/Z);
             cdf2 = (alpha* phi_a - beta* phi_b)/Z;
             %Calculate std deviation of the doubly truncated random variable utilized for rating inference s_{ij}[t],
             %evaluated at x0. 
             f0 = (x0^2*( 1 + cdf2 - cdf1^2));%Eq. 30
             f0_der = 2*x0*( 1 + cdf2 - cdf1^2);
             %------------------------------------------------
             %A better approximation of std. deviation
             y=x0-(f0/f0_der);%Eq. 31
			%Compute the expectted value of s_{ij}[t]
             Ex = m + x0*(normpdf(alpha,0,1) - normpdf(beta,0,1))/Z;%Eq. 32
             x0=real(y);
             
             nu2_next = y;
             if nu2_next>6 || nu2_next<0.9
                 nu2_next = prev_sig^2;
             end
             nu_next = sqrt(nu2_next);
             nu2 = nu2_next;
             if isnan(nu2)
                 nu2 = prev_sig^2;
             end
         end         
         if  pred_mean < part(1) || pred_mean > part(2)
			%Optimal variance for the whole model at t utilized to calculate the posterior distributions of 
			%the latent vectors
             nu3 = abs(sqrt(nu2)+abs(Ex-pred_mean));% Eq. 34
            
             if isnan(nu3)
                 nu3 = nu2;
             end
			%Compute the posterior distributions of the latent vectors 
             post_Sig3 = pinv(invSig1 + (post_mu4*post_mu4' + post_Sig4)/nu3);%posterior cov. matrix of user i
             post_mu3 = post_Sig3*(Ey*post_mu4/nu3 + v1);%posterior mean vector of user i
             post_Sig4 = pinv(invSig2 + (post_mu3*post_mu3' + post_Sig3)/nu3);
             post_mu4 = post_Sig4*(Ey*post_mu3/nu3 + v2);             
             %keep values for the next time (passs them into the data structures lines 183-189)
             next_pred2 = post_mu3'*post_mu4;%Eq. 28
             pred_mean = Ex;
             if next_pred2 >= part(1) && next_pred2 <= part(2)
                post_Sig1 = post_Sig3;
                post_mu1 = post_mu3;
                post_Sig2 = post_Sig4;
                post_mu2 = post_mu4;   
                flag2 = 1;
             else
                flag2 = 0;
             end        
        end
     end  
    if flag2 == 0
        post_Sig1 = pinv(invSig1 + (post_mu2*post_mu2' + post_Sig2)/nu2);
        post_mu1 = post_Sig1*(Ey*post_mu2/nu2 + v1);
        post_Sig2 = pinv(invSig2 + (post_mu1*post_mu1' + post_Sig1)/nu2);
        post_mu2 = post_Sig2*(Ey*post_mu1/nu2 + v2);    
    end
end


obj_u.mu = post_mu1;
obj_u.Sig = post_Sig1;
obj_u.time = event.time;
obj_v.mu = post_mu2;
obj_v.Sig = post_Sig2;
obj_v.time = event.time;


%sigma_stack(n).sigma = nu;
%sigma_stack(n).sigma
end


