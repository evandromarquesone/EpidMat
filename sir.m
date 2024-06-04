% Material Complementar do Livro: 
% Tópicos Matemáticos Aplicados na Modelagem em Epidemiologia - Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

function dy = sir(t,y)
% ============== Declarando o valor de cada Parâmetro ==================
% Esses serão os valores utilizados para computar as coluções numéricas de cada população do modelo 
n=50000; %50 mil indivíduos
beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
nu=1/60; %(a taxa é 1 dividido pelo tempo em que se fica infectado)
%% ================ Equações do Modelo =================================
dy = zeros(3,1); % declarando um vetor nulo, para evitar que o vetor criado na simulação seja não nulo.

dy(1)=-beta*y(1)*y(2);              %Suscetíveis      
dy(2)=beta*y(1)*y(2)-nu*y(2);       %Infectados   
dy(3)=nu*y(2);                      %Recuperados