target_p = 2000;
num_n = 1;
km = .05/num_n;
gamma_m = .001;
gamma_p = .02;

Ca = 1;
Da = .01;
epsilon = .01;
T = 5000;
max_its = 50000000;
CV_squareds = [];
m_burst_average = 20;
ts1 = [];
ms1 = [];
as1 = [];
ts2 = [];
ms2 = [];
as2 = [];
ts3 = [];
ms3 = [];
as3 = [];
for i = 1:3
    if i == 1
       gamma_a = gamma_m*1;
    elseif i == 2
       gamma_a = gamma_m*10;
    else
       gamma_a = gamma_m*45;
    end
    K = 20;
    kp = target_p*gamma_p/m_burst_average;
    pmburst = 1/m_burst_average;
    m = [m_burst_average];
    p = [2000];
    a = [0];
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
            w7 = gamma_a*a;
            W = sum(w1)+w2+w3+w4+w5+w6+w7;
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
         l = length(ps);
         if i == 1
            ts1 = ts;
            ms1 = ms;
            as1 = as;
         elseif i == 2
            ts2 = ts;
            ms2 = ms;
            as2 = as;
         else
            ts3 = ts;
            ms3 = ms;
            as3 = as;
         end
    end
subplot(3,2,1);
plot(ts1,ms1,'r')
hold on;
plot(ts1,20*ones(length(ts1),1),'--k')
set(gca,'Xticklabel',[]) 
ylabel('dilute mRNA','interpreter','latex')
subplot(3,2,2);
plot(ts1,as1,'r')
set(gca,'Xticklabel',[]) 
ylabel('droplet mRNA','interpreter','latex')

subplot(3,2,3); 
plot(ts2,ms2,'b')
hold on;
plot(ts2,20*ones(length(ts2),1),'--k')
set(gca,'Xticklabel',[]) 
ylabel('dilute mRNA','interpreter','latex')

subplot(3,2,4);
plot(ts2,as2,'b')
set(gca,'Xticklabel',[]) 
ylabel('droplet mRNA','interpreter','latex')

subplot(3,2,5); 
plot(ts3,ms3,'g')
hold on;
plot(ts3,20*ones(length(ts3),1),'--k')
xlabel('Time (minutes)','interpreter','latex')
ylabel('dilute mRNA','interpreter','latex')

subplot(3,2,6);
plot(ts3,as3,'g')
xlabel('Time (minutes)','interpreter','latex')
ylabel('droplet mRNA','interpreter','latex')