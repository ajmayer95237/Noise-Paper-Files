target_p = 2000;
km = .05; %Bursting rate

gamma_m = .001; %Cytoplasmic degradation rate
gamma_p = .02; %Protein decay rate
m_burst_average = 20; % Burst intensity
K = 20; %Phase separation threshold
kp = target_p*gamma_p/K; %Translation rate tuned around threshold
Ca = 1; %Association rate
Da = .01; %Dissociaiton rate
gamma_a = .001; %Droplet phase mRNA decay rate
pmburst = 1/m_burst_average;

T = 5000000;
it = 0;
max_its = 50000000;
trials = 10;
CV_squareds = [];
epsilonit = 1;
epsilons = [0.01, 0.05,0.1, 0.5];
for epsilon = epsilons
for i = 1:50
i
CV_k = [];
for k = 1:trials
m = [20];
p = [2000];
a = [50];
t = 0;
data = zeros(max_its,4);
it = 0;
while t < T
    if it > max_its
      break
    end  
    it = it+1;
    data(it,1) = t;
    data(it,2) = p;
    data(it,3) = m;
    data(it,4) = a;
    m_tot = m+a;
    w = 0;
    w1 = km;
    w2 = gamma_m*m;
    w3 = kp*m;
    w4 = gamma_p*p;
    w5 = Ca*(sqrt((m-K)^2+epsilon)+m-K);
    w6 = Da*a;
    w7 = gamma_a*a;
    W = w1+w2+w3+w4+w5+w6+w7;
    r = rand;
    dt = -log(r)/W;
    t = t+dt;
    flag = 0;
    
        w = w+w1;
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
            continue;
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
      if m > 0
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
    a = a-1; %mRNA in droplet phase decays
   
end
 data = data(1:it,:); 
 ts = data(:,1);
 ps = data(:,2);
 ms = data(:,3);
 as = data(:,4);
 l = length(ps);
 i_start = floor(l/100); 
 %CV_k(k) = var(ps(i_start:l))/mean(ps(i_start:l))^2;
  ps_eq = ps(floor(l/100):it);
  ps_eq_squared = ps_eq.^2;
 ts_eq = ts(floor(l/100):it);
 t2 = ts_eq(2:length(ts_eq));
 t1 = ts_eq(1:(length(ts_eq)-1));
 delta_t = t2-t1;
 mean_p_squared = dot(ps_eq_squared(1:(length(ps_eq_squared)-1)),delta_t)/(ts_eq(end)-ts_eq(1));
 mean_p = dot(ps_eq(1:(length(ps_eq)-1)),delta_t)/(ts_eq(end)-ts_eq(1));
 CV_k(k) = (mean_p_squared-mean_p^2)/mean_p^2;

end
 gamma_a = gamma_a+.001;
 CV_squareds(i,epsilonit) = mean(CV_k);
end 
%plot(1:50,CV_squareds,'k')
%hold on
%ylim([0 1.5*10^(-3)])  
epsilonit = epsilonit+1;
end
keyboard
plot(1:50,transpose(CV_squareds),'k')
xlabel('$\gamma_a/\gamma_m$','interpreter','latex')
ylabel('$CV^2$','interpreter','latex')

 