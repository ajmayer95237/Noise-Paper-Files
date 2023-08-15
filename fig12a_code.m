target_p = 2000;
num_n = 1;
km = .05/num_n;
gamma_m = .001;
gamma_p = .02;

Ca = 1; %Association rate
Da = .01; %Dissociaiton rate
gamma_a = .05; %Droplet phase mRNA decay rate
epsilon = .01;

T = 5000000;
%T = 500;
max_its = 20000000;
CV_squareds = [];
m_burst_average = 20;
trials = 10;
%rng(1);
thetas = [];
sigmas = [];
All_CV_squareds = [];
CV_2_mpt = .05;
for q = 1:4
    q
    sigma = .1;
    CV_2_mpt
for i = 1:10   
    %theta = sigma^2/(2*20^2*CV_2_theory);
    theta = sigma^2/(2*20^2*CV_2_mpt);
    mu = m_burst_average;
    kp = target_p*gamma_p/mu;
    pmburst = 1/m_burst_average;
    CV_i = zeros(trials,1);
    m = [m_burst_average];
    p = [2000];
    a = [0];
    t = 0;
    i
    
    for k = 1:trials
        k
        it = 0;
        t = 0;
        m = [m_burst_average];
        p = [2000];
        a = [0];
        mpt = mu;
        data = zeros(max_its,5);
        while t < T
            if it > max_its
              break
            end  
            it = it+1;
            data(it,1) = t;
            data(it,2) = p;
            data(it,3) = m;
            data(it,4) = a;
            data(it,5) = mpt;
            w = 0;
            w1 = km*ones(num_n,1);
            w2 = gamma_m*m;
            w3 = kp*m;
            w4 = gamma_p*p;
            if mpt < 0
               mpt_temp = 0;
            else
               mpt_temp = mpt;
            end
            w5 = Ca*(sqrt((m-mpt_temp)^2+epsilon)+m-mpt_temp);
            w6 = Da*a;
            w7 = gamma_a*a;
            W = sum(w1)+w2+w3+w4+w5+w6+w7;
            r = rand;
            dt = -log(r)/W;
            t = t+dt;
            mpt = mpt+theta*(mu-mpt)*dt+sigma*sqrt(dt)*randn;
            
            flag = 0;
            for j = 1:num_n
                w = w+w1(j);
                if r < w/W %Transcription of mRNA (in a burst)
                    rg = rand;
                    wg = 0;
                    itg = 1;
                    while wg < rg
                          wg = (1 - (1-pmburst)^(itg));
                          itg = itg + 1;
                    end
                    B = itg-1;
                    m = m+B;
                    flag = 1;

                    break
                end 
            end
            if flag == 1
              continue
            end
            w = w+w2;
            if r < w/W %decay of mRNA
              m = m-1;
              continue
            end 
            w = w +w3;  
            if r < w/W %translation of mRNA to protein
               p = p+1;
               continue
            end 
            w = w +w4;
            if r < w/W %decay of protein
              p = p-1;
              continue
            end 
            w = w +w5;
            if r < w/W %mRNA moves to droplet phase
              if m >= 1
                 m = m-1;
                 a = a+1;
              end
              continue
            end 
            w = w+w6;
            if r < w/W
            a = a-1; %mRNA moves to dilute phase
            m = m+1;
            continue
            end
            a = a-1;
         end
         data = data(1:it,:); 
         ts = data(:,1);
         ps = data(:,2);
         ms = data(:,3);
         as = data(:,4);
         mpt = data(:,5);
         l = length(ps);
         ps_eq = ps(floor(l/100):it);
         ps_eq_squared = ps_eq.^2;
         ts_eq = ts(floor(l/100):it);
         t2 = ts_eq(2:length(ts_eq));
         t1 = ts_eq(1:(length(ts_eq)-1));
         delta_t = t2-t1;
         mean_p_squared = dot(ps_eq_squared(1:(length(ps_eq_squared)-1)),delta_t)/(ts_eq(end)-ts_eq(1));
         mean_p = dot(ps_eq(1:(length(ps_eq)-1)),delta_t)/(ts_eq(end)-ts_eq(1));
         CV_i(k) = (mean_p_squared-mean_p^2)/mean_p^2;
    end
    CV_squareds(i) = mean(CV_i);
    sigma = sigma+.1;
end
CV_2_mpt = CV_2_mpt+.05;
All_CV_squareds = [All_CV_squareds; CV_squareds];
end


%%% leg1 = legend('Orstein-Uhlenbeck Phase Separation Threshold Simulation','Deterministic System Driven By Orstein-Uhlenbeck mRNA abundance: $CV^2 = \frac{\gamma_p \sigma^2}{2 \mu^2 \theta (\gamma_p + \theta)}$','Unregulated System: $CV^2 = \frac{\gamma_p \langle B \rangle}{(\gamma_m + \gamma_p) \mu}$')
%%%  set(leg1,'Interpreter','Latex');