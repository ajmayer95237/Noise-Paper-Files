target_p = 2000;
num_n = 1;
km = .05/num_n;
gamma_m = .05;
gamma_p = .02;

Ca = 1;
Da = .01;
epsilon = .01;
T = 10000;
max_its = 50000000;
CV_squareds = [];
m_burst_average = 0;
ts1 = [];
ms1 = [];
ts2 = [];
ms2 = [];
ts3 = [];
ms3 = [];
for i = 1:3
    if i == 1
       m_burst_average = 4;
    elseif i == 2
       m_burst_average = 20;
    else
       m_burst_average = 38;
    end
    K = 20;
    kp = target_p*gamma_p/m_burst_average;
    pmburst = 1/m_burst_average;
    m = [m_burst_average];
    p = [2000];
    a = [200];
    t = 0;
   
        it = 0;
        t = 0;
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
            w5 = Ca*(sqrt((m-K)^2+epsilon)+m-K);
            w6 = Da*a;
            W = sum(w1)+w2+w3+w4+w5+w6;
            r = rand;
            dt = -log(r)/W;
            t = t+dt;
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
            a = a-1; %mRNA moves to dilute phase
            m = m+1; 
         end
         data = data(1:it,:); 
         ts = data(:,1);
         ps = data(:,2);
         ms = data(:,3);
         as = data(:,4);
         l = length(ps);
         if i == 1
            ts1 = ts;
            ms1 = ms;
         elseif i == 2
            ts2 = ts;
            ms2 = ms;
         else
            ts3 = ts;
            ms3 = ms;
         end
    end
subplot(3,1,1);
plot(ts1,ms1,'r')
hold on;
plot(ts1,20*ones(length(ts1),1),'--k')
set(gca,'Xticklabel',[]) 
ylabel('dilute mRNA','interpreter','latex')
subplot(3,1,2); 
plot(ts2,ms2,'b')
hold on;
plot(ts2,20*ones(length(ts2),1),'--k')
set(gca,'Xticklabel',[]) 
ylabel('dilute mRNA','interpreter','latex')

subplot(3,1,3); 
plot(ts3,ms3,'g')
hold on;
plot(ts3,20*ones(length(ts3),1),'--k')
xlabel('Time (minutes)','interpreter','latex')
ylabel('dilute mRNA','interpreter','latex')
