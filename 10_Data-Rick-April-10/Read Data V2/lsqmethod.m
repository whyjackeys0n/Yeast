clear
clc

load S288C_0410_1pct_glucose_pbprch
load S288C_0410_1pct_sucrose_pbprch
load S288C_0410_diff_sucrose_pbpr

unit = 1e6 * 4.6e-11;

pop_density_gL = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_pop_density.*unit;
time = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_time;
conc_pct = S288C_0410_1pct_glucose_pbprch.S288C_0410_1pct_glucose_pbprch_conc_pct;

p1 = [];

figure(1)
for i = 1:length(conc_pct)
    ytrue = pop_density_gL(i,:);
    f = @(c,x)c(1)./(1+exp(-c(2)*(time-c(3))));
    c = lsqcurvefit(f,[5 1 12],time,ytrue);
    p1 = [p1;c];
    if ceil(i/7)==1
        plot(time,ytrue,'o-',time,f(p1(i,:),time),'r.-');
        hold on
    elseif ceil(i/7)==2
        plot(time,ytrue,'o-',time,f(p1(i,:),time),'b.-');
        hold on
    elseif ceil(i/7)==3
        plot(time,ytrue,'o-',time,f(p1(i,:),time),'g.-');
        hold on
    end
end

clear pop_density_gL time conc_pct

pop_density_gL = S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_pop_density.*unit;
time = S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_time;
conc_pct = S288C_0410_1pct_sucrose_pbprch.S288C_0410_1pct_sucrose_pbprch_conc_pct;

p2 = [];

figure(2)
for i = 1:length(conc_pct)
    ytrue = pop_density_gL(i,:);
    f = @(c,x)c(1)./(1+exp(-c(2)*(time-c(3))));
    
    if ceil(i/7)==1
        c = lsqcurvefit(f,[5 1 38],time,ytrue);
        p2 = [p2;c];
        plot(time,ytrue,'o-',time,f(p2(i,:),time),'r.-');
        hold on
    elseif ceil(i/7)==2
        c = lsqcurvefit(f,[5 1 12],time,ytrue);
        p2 = [p2;c];
        plot(time,ytrue,'o-',time,f(p2(i,:),time),'b.-');
        hold on
    elseif ceil(i/7)==3
        c = lsqcurvefit(f,[5 1 12],time,ytrue);
        p2 = [p2;c];
        plot(time,ytrue,'o-',time,f(p2(i,:),time),'g.-');
        hold on
    end
end

mu_private_1 = mean(p1(1:7,2));
mu_cheat_1 = mean(p1(8:14,2));
mu_public_1 = mean(p1(15:21,2));

mu_private_2 = mean(p2(1:7,2));
mu_cheat_2 = mean(p2(8:14,2));
mu_public_2 = mean(p2(15:21,2));

disp(['private growth rate 1: ',num2str(mu_private_1)]);
disp(['cheat growth rate 1: ',num2str(mu_cheat_1)]);
disp(['public growth rate 1: ',num2str(mu_public_1)]);

disp(['private growth rate 2: ',num2str(mu_private_2)]);
disp(['cheat growth rate 2: ',num2str(mu_cheat_2)]);
disp(['public growth rate 2: ',num2str(mu_public_2)]);
