function funPRT_ML(P)

P
r = P.r;
p = P.p;

% loading the details for the forecast evaluation
StartEst = P.StartEst;
StartEv  = P.StartEv;
EndEv    = P.EndEv ;
fcstH_Q = P.fcstH_Q;

% maximum number of EM iterations
P.max_iter = 300;

%--------------------------------------------------------------------------
% Loading monthly data
%--------------------------------------------------------------------------
DataFile = P.DataFile;
LegendFile = P.Legend;

[a,xyz,b] = xlsread(DataFile,'MonthlyLong');

DataM = a(5:end,:);
% if strcmp(version('-release'),'2006b')
DatesM = datestr(cell2mat(b(5:end,1)) + datenum('30-Dec-1899'));  % datenum(cell2mat(b(5:end,1)),'dd/m/yy');
% else
%     DatesM = datenum(b(4:end,1));
% end
DatesMV = datevec(DatesM);
DatesMV = DatesMV(:,1:2);

T = length(DatesM);

[a,xyz,b] = xlsread(LegendFile,'TransfM');
b = b(2:end,:);
GroupM = string(b(:,3));
SeriesM = string(b(:,2));
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
[a,xyz,b] = xlsread(DataFile,'QuarterlyLong');

DataQ = a(5:end,:);

[a,xyz,b] = xlsread(LegendFile,'TransfQ');
b = b(2:end,:);
GroupQ = string(b(:,3));
SeriesQ = string(b(:,2));
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

[a,xyz,b] = xlsread(LegendFile,'Models');
ListIdx = 3; % find(ismember(b(1,3:end),P.Model));
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
% nQ = length(intersect(SeriesQ, Series))
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
% out-of-sample evaluation
%--------------------------------------------------------------------------

iS = find(DatesV(:,1) == StartEv(1) & DatesV(:,2) == StartEv(2));
iE = find(DatesV(:,1) == EndEv(1) & DatesV(:,2) == EndEv(2));

FcstQ = nan(iE,(fcstH_Q+2),nQ);

Month = mod(DatesV(:,2),3);
Month(Month == 0) = 3;

for i = iS:iE

    Date_i = DatesV(i,:);
    Month_i = Month(i);
    disp(['Computing the predictions for the vintages: y', num2str(DatesV(i,1)),' m',num2str(DatesV(i,2))])


    eval(['UnbP = UnbPattM',int2str(Month_i),';']);

    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    n_nan = max([0,(fcstH_Q+1)*3-Month_i+nn]);
    X = [X;nan(n_nan,nVar)];
    % estimation and forecast
    eval(['Res = EM_DFM_SS',P.idio,P.restr,'(X,P);'])
    % collecting the forecast for the quarterly variables
    FcstQ(i,:,:) = Res.X_sm([i-Month_i:3:i-Month_i+3*(fcstH_Q+1)],nM+1:end);

end

FcstQQ = [];
TrueQQ = [];

for k = iS:3:iE-((fcstH_Q+2)*3-1)
    f_t = [];
    for j = 1:fcstH_Q+2
        f_t = cat(2,f_t,FcstQ(k+3*(j-1),fcstH_Q+3-j,:),...
            FcstQ(k+3*(j-1)+1,fcstH_Q+3-j,:),...
            FcstQ(k+3*(j-1)+2,fcstH_Q+3-j,:));
    end
    FcstQQ = cat(1,FcstQQ,f_t);
    TrueQQ = [TrueQQ;Data(k+(fcstH_Q+1)*3-1,nM+1:end)];
end
DateQQ = Dates(iS+(fcstH_Q+1)*3-1:3:iE-3);
DateQQ_V = DatesV(iS+(fcstH_Q+1)*3-1:3:iE-3,:);

% RMSFE for GDP
nF = size(FcstQQ,2);
iVar = find(ismember(SeriesQ,'GDP'));
rmsfeGDP  = sqrt(mean((FcstQQ(1:end,:,iVar) - repmat(TrueQQ(1:end,iVar),1,nF)).^2)');

%--------------------------------------------------------------------------
% saving the results
%--------------------------------------------------------------------------
datafile = ['fcst',P.Model(1:2),P.idio,strrep(int2str(P.r),' ',''),int2str(P.p)];
eval(['save ',datafile,...
    ' FcstQQ TrueQQ DateQQ DateQQ_V Data Series SeriesQ Group Dates DatesV P iS iE rmsfeGDP'])




