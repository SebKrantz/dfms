function funPRT_ML(P)

P
r = P.r;
p = P.p;

% loading the details for the forecast evaluation
StartEst = P.StartEst;
StartEv  = P.StartEv;
EndEv    = P.EndEv ;
Qnews = P.Qnews;
SerNews = P.SerNews;

% maximum number of EM iterations
P.max_iter = 300;

%--------------------------------------------------------------------------
% Loading monthly data
%--------------------------------------------------------------------------
DataFile = P.DataFile;
LegendFile = P.Legend;

[a,b] = xlsread(DataFile,'MonthlyLong');

DataM = a(5:end,:);
% if strcmp(version('-release'),'2006b')
DatesM = datenum(b(4:end,1),'dd/m/yy');
% else
%     DatesM = datenum(b(4:end,1));
% end
DatesMV = datevec(DatesM);
DatesMV = DatesMV(:,1:2);

T = length(DatesM);

[a,b] = xlsread(LegendFile,'TransfM');
b = b(2:end,:);
GroupM = b(:,3);
SeriesM = b(:,2);
%Transformation
TransfM = a(:,6:7);
% unbalancedeness patterns
UnbM = a(:,8:10); 

% Transforming the data
DataMTrf = DataM;
DataMTrf(:,TransfM(:,1) == 1) = 100*log(DataMTrf(:,TransfM(:,1) == 1));
DataMTrf(2:end,TransfM(:,2) == 1) = (DataMTrf(2:end,TransfM(:,2) == 1) ...
    - DataMTrf(1:end-1,TransfM(:,2) == 1));
DataMTrf(1,TransfM(:,2) == 1) = nan;

[tM,nM] = size(DataMTrf);
DataMTrf = [DataMTrf;nan(T-tM,nM)];

%--------------------------------------------------------------------------
% Loading quarterly data
%--------------------------------------------------------------------------
[a,b] = xlsread(DataFile,'QuarterlyLong');

DataQ = a(5:end,:);

[a,b] = xlsread(LegendFile,'TransfQ');
b = b(2:end,:);
GroupQ = b(:,3);
SeriesQ = b(:,2);
%Transformation
TransfQ = a(:,6:7);
% unbalancedeness patterns
UnbQ = a(:,8:10); 

% Transforming the data
DataQTrf = DataQ;
DataQTrf(:,TransfQ(:,1) == 1) = 100*log(DataQTrf(:,TransfQ(:,1) == 1));
DataQTrf(2:end,TransfQ(:,2) == 1) = (DataQTrf(2:end,TransfQ(:,2) == 1) ...
    - DataQTrf(1:end-1,TransfQ(:,2) == 1));
DataQTrf(1,TransfQ(:,2) == 1) = nan;

% quarterly at monthly frequency
DataQMTrf = kron(DataQTrf,[nan;nan;1]);

[tQ,nQ] = size(DataQMTrf);
DataQMTrf = [DataQMTrf;nan(T-tQ,nQ)];

% Building the matrix of restrictions on the loadings for the quarterly
% variables
P.Rconstr = [2 -1 0 0 0;...
    3 0 -1 0 0;...
    2 0 0 -1 0;...
    1 0 0 0 -1];
P.q = zeros(4,1);
P.restr = '_restrMQ';
%--------------------------------------------------------------------------
% complete dataset
%--------------------------------------------------------------------------
Data = [DataMTrf DataQMTrf];
Series = [SeriesM;SeriesQ];
Group = [GroupM;GroupQ];
UnbPatt = [UnbM;UnbQ];

iEst = find(DatesMV(:,1) == StartEst(1) & DatesMV(:,2) == StartEst(2));

[a,b] = xlsread(LegendFile,'Models');
ListIdx = find(ismember(b(1,3:end),P.Model));
List = find(a(:,ListIdx+2)==1);
nM = sum(List<=nM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = Data(iEst:end,List);         %
Dates = DatesM(iEst:end,:);      %
DatesV = DatesMV(iEst:end,:);     %
Series = Series(List);
Group = Group(List);
UnbPatt = UnbPatt(List,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nVar = size(Data,2);
nQ = nVar-nM;

%--------------------------------------------------------------------------
% unbalancedness patterns
%--------------------------------------------------------------------------
nn = min(min(UnbPatt));
nn = min(nn,0);

UnbPattM1 = zeros(12-nn,nVar);
UnbPattM2 = zeros(12-nn,nVar);
UnbPattM3 = zeros(12-nn,nVar);

nUnb = 12-nn;

for i = 1:nVar
    UnbPattM1(end-UnbPatt(i,1)+1+nn:end,i) = nan;
    UnbPattM2(end-UnbPatt(i,2)+1+nn:end,i) = nan;
    UnbPattM3(end-UnbPatt(i,3)+1+nn:end,i) = nan;
end

P.nQ = nQ;

%--------------------------------------------------------------------------
% computing the forecasts and news
% forecast updates are performed each month
%--------------------------------------------------------------------------

iS = find(DatesV(:,1) == StartEv(1) & DatesV(:,2) == StartEv(2));
iE = find(DatesV(:,1) == EndEv(1) & DatesV(:,2) == EndEv(2));
iQ = find(DatesV(:,1) == Qnews(1) & DatesV(:,2) == Qnews(2));
iSer = find(ismember(Series,SerNews));


Month = mod(DatesV(:,2),3);
Month(Month == 0) = 3;

% data at the start of the forecast sequence
Month_i = Month(iS-1);
% first unbalancedeness pattern
eval(['UnbP = UnbPattM',int2str(Month_i),';']);
X = Data(1:iS-1-nn,:);
temp = X(end-nUnb+1:end,:);
temp(isnan(UnbP)) = nan;
X(end-nUnb+1:end,:) = temp;
X_old = [X;nan(max(0,iQ-(iS-1-nn)),nM+nQ)];


DatesNews = Dates(iS:iE);
y_old = zeros(iE-iS+1,1);
y_new = zeros(iE-iS+1,1);
groupnews = zeros(iE-iS+1,length(unique(Group)));
singlenews = zeros(iE-iS+1,nM+nQ);

for i = iS:iE
       
    Date_i = DatesV(i,:);
    Month_i = Month(i);
    disp(['Computing the news for the vintages: y', num2str(DatesV(i,1)),' m',num2str(DatesV(i,2))])
    
    
    % first unbalancedeness pattern
    eval(['UnbP = UnbPattM',int2str(Month_i),';']);
    
    % updated data set
    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    X_new = [X;nan(max(0,iQ-(i-nn)),nM+nQ)];
    
    T_o = size(X_old,1);
    T_n = size(X_new,1);
    X_old = [X_old;nan(T_n-T_o,nM+nQ)];
 
    % new estimates and forecasts
    eval(['R_new = EM_DFM_SS',P.idio,P.restr,'(X_new,P);'])
    R_new.Groups = Group;
    R_new.Series = Series;
    
    % news
    [y_old(i-iS+1,1),y_new(i-iS+1,1),groupnews(i-iS+1,:),singlenews(i-iS+1,:),...
        gainT,serGainT,actual(:,i-iS+1),fcst(:,i-iS+1)] = ...
        News_DFM_ML(X_old,X_new,R_new,iQ,iSer);
    
    X_old = X_new;

  
end

GroupNames = unique(Group)';
trueSer = Data(iQ,iSer);
DateQQ = Dates(iS:iE);
DateQQ_V = DatesV(iS:iE,:);
check = y_new-y_old-sum(groupnews,2)

%--------------------------------------------------------------------------
% saving the results
%--------------------------------------------------------------------------
datafile = ['news',P.Model(1:2),P.idio,strrep(int2str(P.r),' ',''),int2str(P.p)];
eval(['save ',datafile,...
    ' y_old y_new trueSer DateQQ DateQQ_V Data Series GroupNames P groupnews singlenews fcst actual'])




