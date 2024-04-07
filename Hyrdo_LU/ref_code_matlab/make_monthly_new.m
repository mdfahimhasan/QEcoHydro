function [monthly,Z_MONTHLY] =  make_monthly_new(Z,Dates,option)
 %%
% Z      = P_US(:,4:end);
% Dates  = Dates_US;
% option = 2;
%%
 years       = unique(Dates(:,1));
 years_mult  = repmat(years,[12 1]);
 years_mult  = sort(years_mult);
 month_mult  = repmat([1:12]',[length(years) 1]);

%%
Z_MONTHLY    = NaN(1,size(Z,2));

for  i=1:length(years)
%     i=2
    aux1                 = find(Dates(:,1)==years(i));  
    Z_year               = Z(aux1,:); 
    Dates_year           = Dates(aux1,:); 

    for j= 1:12
%     j=1
    aux2                  = find(Dates_year(:,2)==j);   
    
    if isempty(aux2)
    Z_monthly(j,:)        = NaN(1,size(Z,2));
    
    else

    if option == 1
    Z_monthly(j,:)        = mean(Z_year(aux2,:),1);
    end
    
    if option == 2
    Z_monthly(j,:)        = sum(Z_year(aux2,:),1);
    end
   
    end



    end
    Z_MONTHLY  = vertcat(Z_MONTHLY,Z_monthly);
    clear Z_monthly
end

Z_MONTHLY(1,:) = [];
monthly        = [years_mult month_mult Z_MONTHLY];
