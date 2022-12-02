
target_p = 2000;
num_n = 1;
km = .05/num_n;
kp = 2;
gamma_m = .05;
gamma_p = .02;
Da = 1;
%Ca = 20;
mean_as = 0;
m_burst_average = 20;
p_burst_average = 1;
pmburst = 1/m_burst_average;
ppburst = 1/p_burst_average;
T = 500000;

max_its = 20000000;
CV_squareds = [];

mean_ms = m_burst_average*km/gamma_m;
mean_ps = kp/gamma_p*mean_ms;
trials = 5;
for i = 0:50
i
CV_i = [];
mean_as = i;
Ca = Da*mean_as/mean_ms;
for k = 1:trials
    
m = [mean_ms];
p = [mean_ps];
a = [mean_as];
t = 0;
it = 0;
data = zeros(max_its,4);
while t < T
    if it > max_its
      break
    end  
    it = it+1;
    data(it,1) = t;
    data(it,2) = p;
    data(it,3) = m;
    data(it,4) = a;

    w = 0;
    w1 = km*ones(num_n,1);
    w2 = gamma_m*m;
    w3 = kp*m;
    w4 = gamma_p*p;
    w5 = Ca*m;
    w6 = Da*a;
    W = sum(w1)+w2+w3+w4+w5+w6;
    dt = -log(rand)/W;
    t = t+dt;
    r = rand;
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
        rg = rand;
        wg = 0;
        itg = 1;
        while wg < rg
              wg = (1 - (1-ppburst)^(itg));
              itg = itg + 1;
        end
        B = itg-1;
        
        p = p+B;
        continue
    end 
    w = w +w4;
    if r < w/W %decay of protein
      p = p-1;
      continue
    end 
    w = w +w5;
    if r < w/W %mRNA moves to droplet phase
      m = m-1;
      a = a+1;
      continue
    end 
    a = a-1; %mRNA moves to dilute phase
    m = m+1; 
end
 data = data(1:it,:); 
 ts = data(:,1);
 ps = data(:,2);
 ms = data(:,3);
 as = data(:,4);   
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
CV_squareds(i+1) = mean(CV_i);
end
plot(0:1:50,CV_squareds)
hold on
x = 0:1:50;
plot(x,1/target_p+gamma_p*m_burst_average./(gamma_p*x+(gamma_p+gamma_m)*mean_ms),'--')
leg1 = legend('Simulated Data','$\frac{1}{\overline{\langle p \rangle}}+\frac{B}{\overline{\langle a \rangle} + \left(1+\gamma_p/\gamma_m \right)\overline{\langle m \rangle}}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',17);
xlabel('$\overline{\langle a \rangle}$','Interpreter','latex')
ylabel('$CV^2(p)$','Interpreter','latex')






