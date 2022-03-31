step = 0;
for i = 5:1:20
    for j = 1:1:8
        for k = 5:5:100
            for l = 50:50:400
                step = step + 1;
                if step == 2561
                    disp(['Vgmx',num2str(i),' Vomax',num2str(j),' KLa',num2str(k),' Xi',num2str(l)]);
                end
            end
        end
    end
end
