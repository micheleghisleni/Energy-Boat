%% Model from official data

% drag resistance curve presa con:
%   - displacemente 250 kg
%   - xcog = 1.8 m
% in questo script si cercera di trovare la curva che meglio approssima i
% dati forniti

close all
clear 
clc

%% data
speed_kn = [0:5:20]';            % vettore dati velocità in knots
drag_N   = [0 88 212 384 624]';  % vettore dati drag resistance in Newton
N        = length (speed_kn);    % numero di dati


%% polynomial fitting

% usando tutti i dati
max_order = 8;                  % ordine massimo per il fitting (tunable parameter)
th        = zeros(max_order);   % inizializzazione vettore dei parametri

for n_model = 1:max_order

    PHI = zeros(N,n_model+1);                    % generiamo il vettore phi che raccoglie i nostri dati
   
    for ind = n_model:-1:0                       % from the largest
        PHI(:,n_model+1-ind) = speed_kn.^ind;
    end
    
    th(1:n_model+1, n_model) = (PHI'*PHI)\PHI'*drag_N;
    y_ls(:,n_model)          = PHI*th(1:n_model+1,n_model);

end

% tt = 0:0.1:20;
% NN = 20/0.1;
% yy = zeros (6,5);
% for nn = 1:5
% for ii = 1:NN+1
%     for ind = 5:-1:0
%         xx(ind+1, ii) = tt(ii).^ind;
%     end
%     
%     yy(:,ii,nn) = th(1:6,nn)'*xx(:, ii);
%     end
% end
%     


figure (1), plot (speed_kn, drag_N, '.', 'MarkerSize', 20),
hold on
plot (speed_kn, y_ls(:,1:5)); legend ('data', '1 order', '2 order', '3 order', '4 order', '5 order'); xlabel 'Kn', ylabel 'N';
hold off

%% leave one out cross validation

speed_kn_loo = speed_kn;
drag_N_loo   = drag_N;
max_order    = 8;
th_loo       = zeros(max_order,5, max_order);
% Eval(1:max_order).val(1:5) = zeros(8,1);
s_in(:,1)    = zeros (8,1);

for n_model = 1:max_order % questo ciclo si puo fare a matrici
        
    for ii = 1:5
        speed_kn_loo(ii) = [];                          
        drag_N_loo(ii)   = [];
        PHI_loo          = zeros(N-1, n_model+1);

            for ind = n_model:-1:0 % from the largest
                PHI_loo(:,n_model+1-ind) = speed_kn_loo.^ind;
                s_in(n_model-ind+1,1) = speed_kn(ii).^ind;
            end
        
        th_loo(1:n_model+1,ii, n_model) = (PHI_loo'*PHI_loo)\(PHI_loo'*drag_N_loo); % LS formula
        speed_kn_loo                    = speed_kn;
        drag_N_loo                      = drag_N;
        % y_ls_loo(:,ii,n_model) = PHI_loo*th_loo(1:n_model+1,ii,n_model);
        Err.m_ord(n_model).val(ii) = (drag_N(ii)-th_loo(:, ii, n_model)'*s_in)^2;  % calcolo degli errori di validazione per i singoli modelli
    end
    
    Err.m_ord(n_model).cv = sum(Err.m_ord(n_model).val(:))/N; % calcolo dell'errore di cross validation
end
figure (2), plot (1:max_order, [Err.m_ord.cv]);

fprintf ('il miglior modello e'' il polinomio del 3 ordine\n');

%% Polyval: polynomial fitting (per le note, vedi sotto)
min_std_err = Inf;
fitting =zeros(5,1);

for coeff = 1:max_order
    [p,S,mu] = polyfit(speed_kn,drag_N,coeff); 
    [y_fit,std_err] = polyval(p,speed_kn,S,mu);
    fitting = [fitting y_fit];
    if std_err < min_std_err
        min_std_err = std_err;
        best_grade = coeff;
        best_fit = y_fit;
    end
end
fprintf("Il grado con la std. error minore è: %d",best_grade);
fprintf("\nLa std_err vale:");fprintf("%g\n",min_std_err);

figure ('Name','Polyval')
plot (speed_kn, drag_N, '.', 'MarkerSize', 20)
hold on
plot(speed_kn,fitting(:,2:end));
legend ('data', '1 order', '2 order', '3 order', '4 order', '5 order', '6 order', '7 order', '8 order');
%plot(speed_kn,best_fit);
legend('data','best fitting')
xlabel('Kn'),ylabel('N');


% Note:
 
% polyfit ritorna:
% p = i coefficienti per il polinomio p(x) di grado coeff. 
% è il miglior fit per il least-squares error

% S = serve per ottenere una stima dell'errore con la funzione polyval

% mu = un vettore che serve per scalare e centrare la figura. 
% mu(1) è la media, mu(2) è la deviazione standard dei dati 
  
% polyval ritorna:
% il polinomio centrato e stima l'errore standard nel predire 
% un'osservazione futura (std_err)

% quando il grado del fitting diventa maggiore del numero di dati a
% disposizione, la std_err diventa correttamente Infinita
 
% le parti di codice commentate servono per plottare non solo il miglior 
% fitting, ma tutti i fitting calcolati

% Tutta questa parte dovrebbe essere uguale alla prima parte, ma con
% un'altra funzione dei matlab.
%% Fit: esponenziale
options = fitoptions;
options.Normalize = 'on'; % per centrare il plot. ci sono tante altre possibilità

% Fitting esponenziale di primo e di secondo grado.
e1 = fit(speed_kn,drag_N,'exp1',options);% y = a*exp(b*x)
e2 = fit(speed_kn,drag_N,'exp2',options); % y = a*exp(b*x) + c*exp(d*x)
% il grado massimo è il 2

% Plot
figure ('Name','Exponential')
plot (speed_kn, drag_N,'.', 'MarkerSize', 20);
hold on
plot(e1,'k'),plot(e2,'g');
legend('data','1 ordine','2 ordine')
xlabel('Kn'),ylabel('N');

%% Fit : cftoll

cftool("Fitting.sfit")
% sse = sum squared error
% rmse = root mean square error

%% Cross validation
% dataset = [drag_N,speed_kn];
% classificationLearner
% https://it.mathworks.com/help/stats/cvpartition.html
%%
% dataset = [drag_N,speed_kn];
%classificationLearner
%https://it.mathworks.com/help/stats/cvpartition.html
n = length(drag_N); % n° of observations
p = 0.3; % percentage of the data for the testing dataset
partition = cvpartition(n,'Holdout',p);

% Indici della partizione
indx_train = training(partition);
dragTrain = drag_N(indx_train);
speed_knTrain = speed_kn(indx_train);
indx_test = test(partition);
dragTest = drag_N(indx_test);
speed_knTest = speed_kn(indx_test);

c = fit(dragTrain,speed_knTrain,'poly2');
% Estimate loss using cross-validation
err = crossval('mse',dragTrain,drag_N(indx_train))