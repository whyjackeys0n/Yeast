clear
clc

Tn_unit = [];
Xn_unit = [];

for i = 10:1:20
    for j = 1:1:8
        for k = 0.02:0.02:0.2
            for l = 50:50:400
                [Tn,Xn] = Yeast_function(i,j,k,l);
                Tn_unit = [Tn_unit,Tn];
                Xn_unit = [Xn_unit,Xn];
                disp(['Vgmx',num2str(i),' Vomax',num2str(j),' KLa',num2str(k),' Xi',num2str(l)]);
            end
        end
    end
end

save tempdata Tn_unit Xn_unit