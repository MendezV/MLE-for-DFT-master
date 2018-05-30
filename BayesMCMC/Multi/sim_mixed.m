%Let's write a simple markov chain simulation that will take a collection
%of flies (m and f for male and female) and let them move through their
%environment according to some random preferences. 

%PARAMETERS:
if ~already_have_params
N_m = 50;
N_f = 50;
bins_x = 7;
iter_total = 100000;
iter_eq = 100; %Number of iterations before registering the data to occupations so distribution can reach equilibrium.
V_m_wght = 2;
V_f_wght = -2;
m_f_wght = -0.75;
m_m_wght = 0.25;
f_m_wght = -0.75;
f_f_wght = 0.25;
steps = (N_m+N_f)/25;%Number of random steps to make before registering a new value for occ.
watch_sim = 0;
prdc_bndry = 1;
end

%Bins counting down from left to right, top to bottom.
bins = zeros(bins_x,bins_x);
counter = 1;
for i=1:size(bins,1)
    for j=1:size(bins,2)
        bins(i,j) = counter; counter = counter+1;
    end
end

%Find the possible bins flies can move to.
bin_ids = unique(bins);
bin_ids = bin_ids(bin_ids>0);
nbrs = {};
if prdc_bndry    
    for i=1:length(bin_ids)
        [y,x] = find(bins==bin_ids(i));
        if x>1
            left = bins(y,x-1);  else left = bins(y,size(bins,2));
        end
        if x<size(bins,2)
            right = bins(y,x+1);  else right = bins(y,1);
        end
        if y>1
            top = bins(y-1,x);  else top = bins(size(bins,1),x);
        end
        if y<size(bins,1)
            bottom = bins(y+1,x);  else bottom = bins(1,x);
        end
        nbrs{i} = [left,right,top,bottom];
    end
else
    for i=1:length(bin_ids)
        [y,x] = find(bins==bin_ids(i));
        if x>1
            left = bins(y,x-1);  else left = [];
        end
        if x<size(bins,2)
            right = bins(y,x+1);  else right = [];
        end
        if y>1
            top = bins(y-1,x);  else top = [];
        end
        if y<size(bins,1)
            bottom = bins(y+1,x);  else bottom = [];
        end
        nbrs{i} = [left,right,top,bottom];
    end
end

%ASSIGN VEXATIONS:
V_m_values = V_m_wght*[0 1 2 3 4 5 6];
V_m = [];
for i=1:length(V_m_values)
    V_m = [V_m, V_m_values(i)*ones(1,bins_x)];
end

V_f_values = V_f_wght*[0 1 2 3 4 5 6];
V_f = [];
for i=1:length(V_m_values)
    V_f = [V_f, V_f_values(i)*ones(1,bins_x)];
end

%Tracks bin locations of each fly.
N_bins = max(bins(:));

pos_m = round((N_bins-1)*rand(N_m,1))+1;
pos_f = round((N_bins-1)*rand(N_f,1))+1;

occ_m = zeros(iter_total,N_bins);occ_f = zeros(iter_total,N_bins);

occ_m(1,:) = histcounts(pos_m,[1:N_bins+1]);
occ_f(1,:) = histcounts(pos_f,[1:N_bins+1]);

for iter = 2:iter_total
    for st=1:steps
        %Randomly choose a fly
        if rand(1)<=N_m/(N_m+N_f)
            [loc,fly] = datasample(pos_m,1);
            choices = [nbrs{loc} loc];
            delE = zeros(size(choices));
            for i=1:length(delE)
                delE(i) = V_m(choices(i)) - V_m(loc) + m_m_wght*(occ_m(iter-1,choices(i)) - occ_m(iter-1,loc)) + m_f_wght*(occ_f(iter-1,choices(i)) - occ_f(iter-1,loc));
            end
            prob = 1./(1+exp(delE));
            %prob = [1./(1+exp(-delE)), 1/2]; %1 is for staying in the same bin.
            pos_m(fly) = datasample(choices,1,'Weights',prob);
            %prob = prob./(sum(prob));
            %prob = [0, cumsum(prob)];
            %new_bin_ind = find(histcounts(rand(1),prob));
            %bns = [nbrs{loc}, loc];
            %pos_m(fly) = bns(new_bin_ind);
        else
            [loc,fly] = datasample(pos_f,1);
            choices = [nbrs{loc} loc];
            delE = zeros(size(choices));
            for i=1:length(delE)
                delE(i) = V_f(choices(i)) - V_f(loc) + f_m_wght*(occ_m(iter-1,choices(i)) - occ_m(iter-1,loc)) + f_f_wght*(occ_f(iter-1,choices(i)) - occ_f(iter-1,loc));
            end
            prob = 1./(1+exp(delE));
            pos_f(fly) = datasample(choices,1,'Weights',prob);
            %prob = prob./(sum(prob));
            %prob = [0, cumsum(prob)];
            %new_bin_ind = find(histcounts(rand(1),prob));
            %bns = [nbrs{loc}, loc];
            %pos_f(fly) = bns(new_bin_ind);
        end
    end
    occ_m(iter,:) = histcounts(pos_m,[1:N_bins+1]);
    occ_f(iter,:) = histcounts(pos_f,[1:N_bins+1]);
    if watch_sim == 1
        occ_im_m = bins;  occ_im_f = bins;
        for i=1:N_bins
            occ_im_m(occ_im_m==i) = occ_m(iter,i);
            occ_im_f(occ_im_f==i) = occ_f(iter,i);
        end
        subplot(1,2,1); imagesc(occ_im_m);title('Males');axis square
        subplot(1,2,2); imagesc(occ_im_f);title('Females'); axis square
        drawnow
    end
end

occ_m = occ_m([iter_eq:end],:);
occ_f = occ_f([iter_eq:end],:);
        
    
    

        