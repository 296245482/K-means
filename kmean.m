function kmean(choose)

%init some parameters
TOL = 1e-3;
colors = {'r*', 'g*', 'b*', 'y*'};

%load data
switch(choose)
    case 1
        fid = fopen('iris.data');
        trainingData = textscan(fid, '%f%f%f%f%s', 'delimiter', ',');
        x = [trainingData{1,2},trainingData{1,4},trainingData{1,1},trainingData{1,3}];
        K = 3;
        [m, n] = size(x);
        %min-max normalization
        for i=1:n
            x(:,i) = (x(:,i)-min(x(:,i)))/(max(x(:,i))-min(x(:,i))); 
        end
        %show the original data
        figure;
        for i=1:K
            hold on, plot(x((i-1)*50+1:(i-1)*50+50, 1), x((i-1)*50+1:(i-1)*50+50, 2), colors{i});
        end
    case 2
        load('seeds_dataset.txt');
        x = seeds_dataset;
        x(:,8) = [];
        [m, n] = size(x);
        K = 3;
        %min-max normalization
        for i=1:n
            x(:,i) = (x(:,i)-min(x(:,i)))/(max(x(:,i))-min(x(:,i))); 
        end
        x(:,[1;3]) = x(:,[3;1]);
        %show the original data
        figure;
        for i=1:K
            hold on, plot(x((i-1)*70+1:(i-1)*70+70, 1), x((i-1)*70+1:(i-1)*70+70, 2), colors{i});
        end
end


%preprocess
%log normalization
% x(:,1) = log(x(:,1))/log(max(x(:,1)));
% x(:,2) = log(x(:,2))/log(max(x(:,2)));

% give up outliers
for i=1:K
    clusterCenter = zeros(1,n);
    l2distence = zeros(round(m/K),1);
    for j=1:round(m/K)
        clusterCenter = clusterCenter + x((i-1)*50+j,:);
    end
    clusterCenter = clusterCenter/(round(m/K));
    for j=1:round(m/K)
        l2distence(j,1) = norm(x((i-1)*50+j,:)-clusterCenter,2);
    end
    for j=1:round(m/K)
        if l2distence(j,1) > (mean(l2distence)+10*var(l2distence)) || l2distence(j,1) < (mean(l2distence)-10*var(l2distence))
            x((i-1)*50+j,1) = 3;
        end
    end
end
id = x(:,1)==3;
x(id,:)=[];
[m, n] = size(x);

figure;
for i=1:m
    hold on, plot(x(i, 1), x(i, 2), 'b+');
end

label = ones(m,1);

%choose K point to begin
% rand = randperm(round(m/K));
% for i=1:K
%     center(i,:) = x(rand(1)+(i-1)*round(m/K),:);
% end
rand = randperm(m);
center=x(rand(1:K),:);
firstCenter = center;

oldSSE = 0;
count = 1;
while 1
    %divide point into K closest
    for i=1:m
        minIndex = 1;
        minDis = norm(x(i,:)-center(minIndex,:),2);
        for j=1:K
            newDis = norm(x(i,:)-center(j,:),2);
            if newDis < minDis
                minDis = newDis;
                minIndex = j;
            end
        end
        label(i) = minIndex;
    end
    
    %find new center
    for i=1:K
        center(i,:) = sum(x(find(label==i),:));
        center(i,:) = center(i,:) / length(find(label==i));
    end
    
    %calculate SSE
    newSSE = 0;
    for i=1:m
        newSSE = newSSE + norm(x(i,:)-center(label(i),:),2);
    end
    fprintf('==========iter %d==========\n',count);
    fprintf('newSSE is %d, SSE is %d\n',newSSE,oldSSE);
    fprintf('improve is %d\n\n',abs(oldSSE-newSSE));
    
    %compara two SSE
    if abs(oldSSE-newSSE) < TOL
        break;
    end
    oldSSE = newSSE;
    count = count + 1;
end

% show plot of clustering
figure;
for i=1:K
    hold on, plot(x(find(label == i), 1), x(find(label == i), 2), colors{i});
end

hold on;
plot(firstCenter(1:K,1),firstCenter(1:K,2),'y*');

end