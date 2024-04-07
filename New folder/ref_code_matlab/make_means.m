function [INPUT_MEAN_MONTHLY, INPUT_MEAN_ANNUAL] =  make_means(INPUT,Dates)


years = unique(Dates(:,1));
years(1) = [];
for i =1:length(years)
        aux                  = find(Dates(:,1) == years(i));
        INPUT_ANNUAL(i,:)    = nansum(INPUT(aux,:));
       
end
INPUT_MEAN_ANNUAL   = nanmean(INPUT_ANNUAL,1);


for i =1:12
        aux                       = find(Dates(:,2) == i);
        INPUT_MEAN_MONTHLY(i,:)   = nanmean(INPUT(aux,:));
       
end
