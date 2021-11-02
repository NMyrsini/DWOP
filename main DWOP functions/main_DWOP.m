%Dynamic weighted ordered probit model integrated to Collaborative Kalman filtering for hotel rating prediction
clc
clear memory
clear all
%tic
showInfo=true;
skip = true;

if (skip)
    %Load Data to construct structure
    if (showInfo)
        disp('Reading files started...')
    end
    %tripadvisor_dataset
    %!! Please insert the corret file path
    fid1 = fopen('tripadvisor_dataset\Preveza_Tripadvisor.txt');
    C = textscan(fid1, '%d %d %f %s ', 'delimiter',',');
    fclose(fid1);
    if (showInfo)
       disp('finished.')
    end
    elements=numel(C{1});
	%5 segments of the ordered probit model, one for each star (5-star rating system)
    num_part = 5;

    times = single(datenum(C{4},'yyyy-mm-dd'));
    Means = single(zeros(elements,1));

    %__________Hash map construction (instead of pointers)__________________
    revMapUsers = containers.Map;
    users = uint32(zeros(elements,1));
    userCount = single(1);

    if (showInfo)
        disp('Events structure creation Started...')
    end

    numOfHotels1 =numel(unique(C{1}));
    numOfHotels =max(C{1});
    %Latent state vector dimensionality 
    d=5;a=1;b=5;
	%Initialize the covariance matrix with ones on the main diagonal and zeros elsewhere
    SigVal = single(eye(d,d));

    counter = 0;counter2 = 0;
    if (showInfo)
        disp('finished.')
        %toc
    end
    % User & hotel structure declaration
    numOfUsers = numel(unique(C{2}));
    %User latent vector
    obj_u = repmat(struct('mu',0,'Sig',SigVal,'time',single(datenum(C{4}(1),'yyyy-mm-dd')),'alpha',0,'beta',1),1,numOfUsers);
    %Hotel latent vector
	obj_v = repmat(struct('mu',0,'Sig',SigVal,'time',single(datenum(C{4}(1),'yyyy-mm-dd')),'beta',0.1),1, numOfHotels);
    % obj_u (user) & obj_v (hotel) mean vectors random initialization
    for x =1:numel(obj_u)
       obj_u(x).mu = single(sqrt(a/d) + sqrt((b-a)/d).*rand(d,1));
    end

    for jj = 1:numel(obj_v)
       obj_v(jj).mu = single(sqrt(a/d) + sqrt((b-a)/d).*rand(d,1));
    end
	%Every Event represents a time stamp in the time series of ratings which embeds user id, hotel id, time step, etc.
    events = repmat(struct('id1', 0, 'id2',0, 'z',0,'time',0,'mean',0,'var',0,'lbound',0,'rating',3),1, elements);
    %Structure of events initialization
    for i =1:elements
        events(i).id1=uint32(0);
        events(i).id2=uint16(0);
        events(i).z=uint8(0);
        events(i).time=single(0);
        events(i).mean=single(0);
        events(i).var=single(0);
        events(i).lbound=single(0);
    end

    %___________Hash-map instead of pointers_____________
    for i=1:elements
            counter = counter+1;
            %___Hash map______________________________________________________
            if ~isKey(revMapUsers,num2str(C{2}(counter)))
                revMapUsers(num2str(C{2}(counter))) = userCount; %revMapUsers(num2str(1366860)) = 2 
                users(userCount) = uint32(C{2}(counter));%real user id in ith event
                userCount = userCount + 1;
            end
            events(i).id1 =uint32(revMapUsers(num2str(C{2}(counter))));
            events(i).id2 = uint16(C{1}(counter));
            events(i).z = uint8(C{3}(counter));
            events(i).time = times(counter);

    end
    hot = unique(C{1});
    for kk=1:numel(unique(C{1}))
        ID{kk}=find((C{1}==hot(kk)));
        FirstDates{kk}= ID{kk}(1);
    end
    %initialization to compute beta param. of N(beta,0.1)
    z=zeros(1,numOfHotels);
    summ=zeros(1,numOfHotels);
    summs=zeros(1,numOfHotels);
    mu=zeros(1,numOfHotels);
    beta_x=zeros(1,numOfHotels);

    if (showInfo)
        disp('finished.')
    end

    if (showInfo)
        disp('DWOP started...')
    end
    
    %Start sliding the window over the time series of ratings and apply DWOP until the last rating is met.
    [predictions,obj_u,obj_v,events,counter2, elapsed_time] = sliding_windows(events,obj_u,obj_v,num_part,z,summ,mu,numOfHotels,FirstDates);
    time_avg = mean(elapsed_time)
    time_sum = sum(elapsed_time);
    
    
    if (showInfo)
        disp('finished.')
    end
       %__________Save predicted rating
    for p=1:numel(events)
        Means(p)=events(p).rating;
        %Means(p)=events(p).mean;
    end

    s.(sprintf('PredMean%d', 1)) = single(Means); % Save it to a new variable in a structure
	%seve predictions in a file
    save('evaluation_DWOP\Preveza_Tripadvisor_DWOP_predictions.mat', '-struct', 's'); 
    clear s 
end