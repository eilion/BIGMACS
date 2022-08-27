function [Info] = getInfo(Samples)

L = length(Samples);

Info = struct('name',cell(L,1),'depth',cell(L,1),'lower',cell(L,1),'median',cell(L,1),'upper',cell(L,1));

for ll = 1:L
    Info(ll).name = Samples(ll).name;
    Info(ll).depth = Samples(ll).depth;
    
    path = ['Cores/',Info(ll).name,'/query_depths.txt'];
    
    if exist(path,'file') == 2
        fileID = fopen(path);
        INFO = textscan(fileID,'%s %s');
        fclose(fileID);
        if sum(strcmp(INFO{1},'start_depth:')==1) == 1
            st = str2double(INFO{2}{strcmp(INFO{1},'start_depth:')==1});
            st = max(st,Info(ll).depth(1));
        else
            st = Info(ll).depth(1);
        end
        if sum(strcmp(INFO{1},'end_depth:')==1) == 1
            ed = str2double(INFO{2}{strcmp(INFO{1},'end_depth:')==1});
            ed = min(ed,Info(ll).depth(end));
        else
            ed = Info(ll).depth(end);
        end
        if sum(strcmp(INFO{1},'interval:')==1) == 1
            itvr = str2double(INFO{2}{strcmp(INFO{1},'interval:')==1});
        else
            itvr = 0.01;
        end
        
        Depth = (st:itvr:ed)';
        if abs(Depth-ed) > 1e-15
            Depth = [Depth;ed];
        end
        AA = zeros(length(Depth),size(Samples(ll).ages,2));
        
        for m = 1:size(AA,2)
            AA(:,m) = interp1(Samples(ll).depth,Samples(ll).ages(:,m),Depth);
        end
        
        N = size(AA,1);
        LMU = zeros(N,6);
        
        for n = 1:N
            LMU(n,1) = quantile(AA(n,:),0.025);
            LMU(n,2) = quantile(AA(n,:),0.16);
            LMU(n,3) = quantile(AA(n,:),0.5);
            LMU(n,4) = mean(AA(n,:));
            LMU(n,5) = quantile(AA(n,:),0.84);
            LMU(n,6) = quantile(AA(n,:),0.975);
        end
        
        Info(ll).lower = LMU(:,1:2);
        Info(ll).median = LMU(:,3:4);
        Info(ll).upper = LMU(:,5:6);
        Info(ll).depth = Depth;
    else
        AA = Samples(ll).ages;
        N = size(AA,1);
        LMU = zeros(N,6);
        
        for n = 1:N
            LMU(n,1) = quantile(AA(n,:),0.025);
            LMU(n,2) = quantile(AA(n,:),0.16);
            LMU(n,3) = quantile(AA(n,:),0.5);
            LMU(n,4) = mean(AA(n,:));
            LMU(n,5) = quantile(AA(n,:),0.84);
            LMU(n,6) = quantile(AA(n,:),0.975);
        end
        
        Info(ll).lower = LMU(:,1:2);
        Info(ll).median = LMU(:,3:4);
        Info(ll).upper = LMU(:,5:6);
    end 
end


end