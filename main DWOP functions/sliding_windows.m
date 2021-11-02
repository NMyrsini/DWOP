function [predictions,obj_u,obj_v,events,counter2, elapsed_time] = sliding_windows(events,obj_u,obj_v,num_part,z,summ,mu,numOfHotels,FirstDates)
%
% For the nth event:
%   events(n).id1  : index of objects struct of participant in obj_u
%   events(n).id2  : index of objects struct of participant in obj_v
%   events(n).z    : response in {1,...,num_part}
%   events(n).time : time stamp
%
% For the ith object in group u (similar for group v):
%   obj_u(i).mu   : posterior mean of q distribution
%   obj_u(i).Sig  : posterior covariance of q distribution
%   obj_u(i).time : time stamp of last event


num_events = length(events);
predictions = zeros(1,num_events);
%initialize segments (part)
part = zeros(1,num_part+1);
%prefixed standard deviation
nu = 1.76;
target_sigma = nu;
%initialize the left boundary of the segment that corresponds to 1 star (its value does not play a critical role)
left = -10.8;
%Create a data structure to keep the computed values of std. deviation (sigma)
sigma_stack = repmat(struct('HotelName', 1, 'sigma',0),1, numel(events));
%count execution time (for every window and then the average value will be
%computed in main function)
elapsed_time = zeros(num_events,1);

counter2 =0;
for n = 1:numel(events)%every n corresponds to a time step
    t=tic;
         
        for i=1:numOfHotels 
            if events(n).id2==i
                z(i)=z(i)+1;
                summ(i)=summ(i)+events(n).z;
                m_summ(i) = summ(i)/z(i);
                obj_v(events(n).id2).beta = m_summ(i);
            end
        end
        %sliding window size is 4 
        len = 4;
        if mod(n,len)==0 
			%initialize to 3 stars to confront cold start problem
            obj_v(events(n).id2).beta = 3;
            for k = 1:numOfHotels
                if events(n).id2==k
                    z(k)=0;
                    summ(k)=0;
                    mu(k)=0;
                    m_summ(k) = 0;
                end
            end
        end
      
          
        %keep the average ratings of a hotel (at the previous time steps)
        nus = obj_v(events(n).id2).beta;
        
        
        
        sigma_stack(n).HotelName = events(n).id1; 
        for  jj=1:numel(FirstDates)
            if n == 1 || n==FirstDates{jj}
                initBool = 1;
                sigma_stack(n).sigma = target_sigma;
                break
            else 
                initBool = 0;  
            end
        end 
		%If it is the first a hotel appears set the previous value of the std. deviation equal
%		to the prefixed one to confront the lack of previous information
        if initBool == 1 
            prev_sig =  target_sigma;
            counter2 = counter2+1;
        else
            for j = 2 : n-1 

                if sigma_stack(n).HotelName == sigma_stack(n-j).HotelName
                    prev_sig = sigma_stack(n-j).sigma; 
                    break 
                else
                    prev_sig = target_sigma;
                end  
            end

        end
		nu_=prev_sig*len;%auxiliary buffer
		% nu: std. deviation computed through NR update etc., defining the width of segments at each time step.
        nu = prev_sig;
		%Construct the segments of the ordered probit model. Check if the average rating of a hotel (in the past) falls into a 
		%specific region corresponding to the number of rating stars. The mapping is conducted directly on the desired segment.
        if nus > 3 && nus<= 4
            part(1) = left;
            part(2) = part(1) + nu;
            part(3) = part(2)+ nu;
            part(4) = part(3) + nu;
            part(5) = part(4) + nu;
            part(6) = part(5) + nu;
            nu = nu_;
        elseif nus >4  
            part(1) = left;
            part(2) = part(1) + nu;
            part(3) = part(2) + nu;
            part(4) = part(3) + nu;
            part(5) = part(4) + nu;   
            part(6) = part(5) + nu;
            nu = nu_;
         elseif nus>1 && nus<=2          
            part(1) = left;
            part(2) = part(1) + nu;
            part(3) = part(2) + nu;
            part(4) = part(3) + nu;
            part(5) = part(4) + nu;   
            part(6) = part(5) + nu;
            nu = nu_;
        elseif nus > 2 && nus<=3         
            part(1) = left;
            part(2) = part(1) + nu;
            part(3) = part(2) + nu;
            part(4) = part(3) + nu;
            part(5) = part(4) + nu;   
            part(6) = part(5) + nu;
            nu = nu_;            
        elseif  nus< 2        
            part(1) = left;
            part(2) = part(1) + nu;
            part(3) = part(2) + nu;
            part(4) = part(3) + nu;
            part(5) = part(4) + nu;   
            part(6) = part(5) + nu;
            nu = nu_;  
        end

        %A function which applies NR update to compute the std. deviation, accept-reject sampling to compute the expected value
		%of an auxiliary truncated normal random variable, and CKF to compute the posterior moments of latent vectors, and to predict a rating at the current time step.
        [obj_u(events(n).id1),obj_v(events(n).id2),m,v, nu_next,sigma_stack] = DWOP_update(events(n),obj_u(events(n).id1),obj_v(events(n).id2),[part(events(n).z),part(events(n).z+1)],nu,1,sigma_stack,n,initBool);
        %keep the currect sigma (optimal std. deviation)
        sigma_stack(n).sigma = nu_next;
        %rating prediction
        events(n).mean = m;
        
        %scale back to 5 star sating system to save predictions in an understandable format 
        if events(n).mean >= part(1) && events(n).mean<part(2)
            events(n).rating =1;
        elseif events(n).mean >= part(2) && events(n).mean<part(3)
            events(n).rating =2;
        elseif events(n).mean >= part(3) && events(n).mean<part(4)
            events(n).rating =3;
        elseif events(n).mean >= part(4) && events(n).mean<part(5)
            events(n).rating =4;
        elseif events(n).mean >= part(5) 
            events(n).rating =5; 
        end
        

       
t=toc(t);
elapsed_time(n) = t ;
        
       
end
