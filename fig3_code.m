ts_lowCa = [];
ms_lowCa = [];
ps_lowCa = [];
as_lowCa = [];
ts_highCa = [];
ms_lowCa = [];
ps_lowCa = [];
as_lowCa = [];
for i = 1:2
target_p = 2000;
km = .05;
gamma_m = .05;
gamma_p = .02;
Da = 1;
if i == 1
   Ca = .1;
else
   Ca = 10;
end
m_burst_average = 20;
mean_ms = m_burst_average*km/gamma_m;
kp = gamma_p*target_p/(mean_ms);
pmburst = 1/m_burst_average;
T = 10000;
max_its = 20000000;

m = [mean_ms];
p = [target_p];
a = [0];
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
    w1 = km;
    w2 = gamma_m*m;
    w3 = kp*m;
    w4 = gamma_p*p;
    w5 = Ca*m;
    w6 = Da*a;
    W = w1+w2+w3+w4+w5+w6;
    dt = -log(rand)/W;
    t = t+dt;
    r = rand;

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
 if i == 1
    ts_lowCa = ts;
    ms_lowCa = ms;
    ps_lowCa = ps;
    as_lowCa = as;
 else
    ts_highCa = ts;
    ms_highCa = ms;
    ps_highCa = ps;
    as_highCa = as;
 end
end
subplot(3,1,1);
plot(ts_lowCa,ms_lowCa)
hold on
plot(ts_highCa,ms_highCa)
legend('C_a = 0.1','C_a = 20')
title('Active mRNA')
% Plot 2
subplot(3,1,2);
plot(ts_lowCa,as_lowCa)
hold on
plot(ts_highCa,as_highCa)
legend('C_a = 0.1','C_a = 20')
title('Inactive mRNA')

% Plot 3
subplot(3,1,3);
plot(ts_lowCa,ps_lowCa)
hold on
plot(ts_highCa,ps_highCa)
legend('C_a = 0.1','C_a = 20')
title('Protein')









