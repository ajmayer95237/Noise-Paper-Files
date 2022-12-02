target_p = 2000;
num_n = 1;
km = .05/num_n;
gamma_m = .05;
gamma_p = .02;

Ca = 1;
Da = .01;
epsilon = .01;
T = 5000000;
max_its = 50000000;
CV_squareds = [];
m_burst_average = 0;
trials = 5;
for i = 1:40
    i
    m_burst_average = m_burst_average+1;
    K = 20;
    kp = target_p*gamma_p/m_burst_average;
    pmburst = 1/m_burst_average;
    CV_i = zeros(trials,1);
    m = [m_burst_average];
    p = [2000];
    a = [0];
    t = 0;
    for k = 1:trials
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
    
end
figure(1)
plot(1:40,CV_squareds)

% function createfigure(X1, Y1, X2, Y2)
% %CREATEFIGURE(X1, Y1, X2, Y2)
% %  X1:  vector of x data
% %  Y1:  vector of y data
% %  X2:  vector of x data
% %  Y2:  vector of y data
% 
% %  Auto-generated by MATLAB on 01-Sep-2021 09:27:51
% 
% % Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% % Create plot
% plot(X1,Y1,'LineWidth',1,'Color',[0 0 0]);
% 
% % Create plot
% plot(X2,Y2);
% 
% % Create ylabel
% ylabel('$CV^2$','Interpreter','latex');
% 
% % Create xlabel
% xlabel({'Mean Burst Size'},'Interpreter','latex');
% 
% box(axes1,'on');
% % Create line
% annotation(figure1,'line',[0.519642857142857 0.517857142857143],...
%     [0.106142857142857 0.926190476190476],'LineWidth',1,'LineStyle',':');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.197428571428571 0.233333333333333 0.220428571428572 0.0880952380952399],...
%     'String',{'Poisson limit','$1/\langle p \rangle$ = 0.005'},...
%     'Interpreter','latex',...
%     'FitBoxToText','off');
% 
% % Create arrow
% annotation(figure1,'arrow',[0.2875 0.344642857142857],...
%     [0.229952380952381 0.116666666666667]);
% 
% % Create textarrow
% annotation(figure1,'textarrow',[0.569642857142857 0.521428571428571],...
%     [0.60852380952381 0.564285714285715]);
% 
% % Create ellipse
% annotation(figure1,'ellipse',...
%     [0.191142857142857 0.868105515587523 0.0248571428571429 0.0341726618704996],...
%     'Color',[1 0 0],...
%     'LineWidth',1,...
%     'FaceColor',[1 0 0]);
% 
% % Create ellipse
% annotation(figure1,'ellipse',...
%     [0.506342857142859 0.347721822541962 0.0248571428571429 0.0341726618704997],...
%     'Color',[0 0 1],...
%     'LineWidth',1,...
%     'FaceColor',[0 0 1]);
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.570642857142856 0.563549160671463 0.195428571428571 0.13645083932854],...
%     'String',{'Phase Separation','Threshold','$m_{pt} = 20$'},...
%     'Interpreter','latex',...
%     'FitBoxToText','off');
% 
% % Create ellipse
% annotation(figure1,'ellipse',...
%     [0.851942857142859 0.206235011990405 0.0248571428571428 0.0341726618704997],...
%     'Color',[0 1 0],...
%     'LineWidth',1,...
%     'FaceColor',[0 1 0]);
% 
% 
% 
