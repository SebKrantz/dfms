function [X,indNaN] = remNaNs(X,options);

% procedure to replace missing observations
% used for ontaining the initial parameters for EM iterations for 
% “MAXIMUM LIKELIHOOD ESTIMATION OF FACTOR MODELS ON DATA SETS WITH ARBITRARY PATTERN OF MISSING DATA.”
% by Marta Ba?bura and Michele Modugno 

[T,N]=size(X);
k=options.k;
indNaN=isnan(X);
switch options.method
    case 1  % replace all the missing values by moving averages (with median)
        for i = 1:N  
            x = X(:,i);
            x(indNaN(:,i))=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(indNaN(:,i))=x_MA(indNaN(:,i));
            X(:,i)=x;
        end
    case 2 %replace missing values after removing leading and closing nans 
        rem1=(sum(indNaN,2)==N);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            x(indNaN(:,i))=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(indNaN(:,i))=x_MA(indNaN(:,i));
            X(:,i)=x;
        end
    case 3 %only remove rows with leading and closing nans
        rem1=(sum(indNaN,2)==N);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN = isnan(X);
    case 4 %replace missing values with splines after removing leading and closing nans 
        rem1=(sum(indNaN,2)==N);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            t1 = min(find(~isnanx));
            t2 = max(find(~isnanx));
            x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
            isnanx = isnan(x);
            x(isnanx)=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(isnanx)=x_MA(isnanx);
            X(:,i)=x;
        end
    case 5 %replace missing values with splines
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            t1 = min(find(~isnanx));
            t2 = max(find(~isnanx));
            x(t1:t2) = spline(find(~isnanx),x(~isnanx),(t1:t2)');
            isnanx = isnan(x);
            x(isnanx)=nanmedian(x);
            x_MA =filter (ones(2*k+1,1)/(2*k+1),1,[x(1)*ones(k,1);x;x(end)*ones(k,1)]);
            x_MA=x_MA(2*k+1:end);
            x(isnanx)=x_MA(isnanx);
            X(:,i)=x;
        end
    case 6 %replace missing values with random draws after removing leading and closing nans 
        rem1=(sum(indNaN,2)==N);
        nanLead=(cumsum(rem1)==(1:T)');
        nanEnd=(cumsum(rem1(end:-1:1))==(1:T)');
        nanEnd=nanEnd(end:-1:1);
        nanLE=(nanLead | nanEnd);
        X(nanLE,:)=[];
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            x(isnanx) = randn(sum(isnanx),1);
            X(:,i)=x;
        end
    case 7 %replace missing values with random draws
        indNaN=isnan(X);
        for i = 1:N  
            x = X(:,i);
            isnanx = isnan(x);
            x(isnanx) = randn(sum(isnanx),1);
            X(:,i)=x;
        end
end
