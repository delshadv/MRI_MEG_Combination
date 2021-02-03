function [f1] = plot_results(titles,Acc,pos_titles,c,Acc2)

y_label = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'};
N = size(Acc,1);

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
        ylabel(sprintf([y_label{k} num2str(1)]),'Rotation',0,...
            'fontweight','bold','fontsize',16);
        ytickangle(45);
        yticks([0 size(Acc,1)/5])
        yticklabels({'0',sprintf(num2str(size(Acc,1)/5))})
        title(titles{k});
        axi = gca;
        axi.FontSize = 12;
        axi.Position = axi.Position .* [1 1 0.5 0.8];
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
        ax{k} = subplot(size(a,2),2,2*k-1);
        n = histcounts(a(:,k));
        histogram(a(:,k));
        ylabel(sprintf([y_label{k} num2str(1)]),'Rotation',0,...
            'fontweight','bold','fontsize',16);
        ytickangle(45);
        yticks([0 size(Acc,1)/5])
        yticklabels({'0',sprintf(num2str(size(Acc,1)/5))})
        %       axis('tight');
        title(titles{k});
        axi = gca;
        axi.FontSize = 12;
        
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
    
    %f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    ax1 = cell(1,size(ca,2)); % preallocate
    for k = 1 : size(ca,2)
        ax1{k} = subplot(size(ca,2),2,2*k);%,'position',ax{k}.Position);
        n = histcounts(a(:,k));
        histogram(ca(:,k),'FaceColor',[0.4660 0.6740 0.1880]);
        ylabel(sprintf([y_label{k} num2str(2)]),'Rotation',0,...
            'fontweight','bold','fontsize',16);
        ytickangle(45);
        yticks([0 size(Acc,1)/5])
        yticklabels({'0',sprintf(num2str(size(Acc,1)/5))})
        title(pos_titles{k})
        axi = gca;
        axi.FontSize = 12;
        max_y = ylim;
        line([0 0]',[0 max(n)]','Color','k','LineWidth',2)
    end
    
    linkaxes([ax1{:}],'x')
    linkaxes([ax1{:}],'y')


elseif nargin==5
    
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
        ax{k} = subplot(size(a,2),2,2*k-1);
        n = histcounts(a(:,k));
        histogram(a(:,k));
        ylabel(sprintf([y_label{k} num2str(1)]),'Rotation',0,...
            'fontweight','bold','fontsize',16);
        ytickangle(45);
        yticks([0 size(Acc,1)/5])
        yticklabels({'0',sprintf(num2str(size(Acc,1)/5))})
        %       axis('tight');
        title(titles{k});
        axi = gca;
        axi.FontSize = 12;
        
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
    a = zeros(size(Acc2,1),numel(titles)); % preallocate
        for k = 1 : size(Acc2,2)
            a(:,k) = mean(Acc2(:,k,:),3);
            
        end
        
    ca=a*c';
    
    %CI = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
    %ci= CI(ca,95)
    
    p=[];
    for k = 1:size(ca,2)
        p(k) = length(find(ca(:,k) <= 0))/size(ca,1);
    end
    p
    
    %f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    ax1 = cell(1,size(ca,2)); % preallocate
    for k = 1 : size(ca,2)
        ax1{k} = subplot(size(ca,2),2,2*k);%,'position',ax{k}.Position);
        n = histcounts(a(:,k));
        histogram(ca(:,k),'FaceColor',[0.4660 0.6740 0.1880]);
        ylabel(sprintf([y_label{k} num2str(2)]),'Rotation',0,...
            'fontweight','bold','fontsize',16);
        ytickangle(45);
        yticks([0 size(Acc,1)/5])
        yticklabels({'0',sprintf(num2str(size(Acc,1)/5))})
        title(pos_titles{k})
        axi = gca;
        axi.FontSize = 12;
        %max_y = ylim;
        line([0 0]',[0 max(n)]','Color','k','LineWidth',2)
    end
    
    linkaxes([ax1{:}],'x')
    linkaxes([ax1{:}],'y')
    else
        error("Please check inputs")
end


end

