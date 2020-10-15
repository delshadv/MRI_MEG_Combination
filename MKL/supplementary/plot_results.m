function [f1,f2] = plot_results(titles,Acc,pos_titles,c)

N = size(Acc,1)

if nargin==2
    if numel(titles) == size(Acc,2)
        a = zeros(size(Acc,1),numel(titles)); % preallocate
        for k = 1 : size(Acc,2)
            a(:,k) = mean(Acc(:,k,:),3);
            
        end
    else
        error("Please correct titles")
    end
    
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    ax = cell(1,size(a,2)); % preallocate
    for k = 1 : size(a,2)
        ax{k} = subplot(size(a,2),1,k);
        n = histcounts(a(:,k));
        histogram(a(:,k));
        title(titles{k})
        line([50 50]',[0 max(n)]','Color','k','LineWidth',2)
    end
    
    p=[];
    for k = 1:size(a,2)
        p(k) = length(find(a(:,k) <= 50))/size(a,1);
    end
    p
    
    linkaxes([ax{:}],'x')
    linkaxes([ax{:}],'y')
    
elseif nargin==4
    if numel(titles) == size(Acc,2)
        a = zeros(size(Acc,1),numel(titles)); % preallocate
        for k = 1 : size(Acc,2)
            a(:,k) = mean(Acc(:,k,:),3);
            
        end
    else
        error("Please correct titles")
    end
    
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    ax = cell(1,size(a,2)); % preallocate
    for k = 1 : size(a,2)
        ax{k} = subplot(size(a,2),1,k);
        n = histcounts(a(:,k));
        histogram(a(:,k));
        title(titles{k})
        line([50 50]',[0 max(n)]','Color','k','LineWidth',2)
    end
    
    p=[];
    for k = 1:size(a,2)
        p(k) = length(find(a(:,k) <= 50))/size(a,1);
    end
    p
 
    linkaxes([ax{:}],'x')
    linkaxes([ax{:}],'y')
    
    
    
    %% now any contrasts...
    ca=a*c';
    
    %CI = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
    %ci= CI(ca,95)
    
    p=[];
    for k = 1:size(ca,2)
        p(k) = length(find(ca(:,k) <= 0))/size(ca,1);
    end
    p
    
    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    ax = cell(1,size(ca,2)); % preallocate
    for k = 1 : size(ca,2)
        ax{k} = subplot(size(ca,2),1,k);
        n = histcounts(a(:,k));
        histogram(ca(:,k));
        title(pos_titles{k})
        line([0 0]',[0 max(n)]','Color','k','LineWidth',2)
    end
    
    linkaxes([ax{:}],'x')
    linkaxes([ax{:}],'y')
else
    error("Please check inputs")
end

end

