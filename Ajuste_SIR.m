% Material Complementar do Livro: 
% Introdução à Epidemiologia Matemática: Métodos em Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

%% = Ajuste de Dados de Infectados à Solução de Infectados do Modelo SIR =
% Nesse programa fazemos o uso da função fminsearch, para minimizar a 
% distância entre o valor ajustado e o valor observado. O programa foi
% feito para determinar os valores dos parâmetros de infecção e de 
% recuperação, que otimizam essa distância.  
% Confira a versão atualizada em https://github.com/evandromarquesone/EpidMat

clc;  clear all;  close all; 

% % Todos os dados de casos de indivíduos infectados
 B=[1 3;2 8;3 28;4 75;5 221;6 291;7 255;8 235;9 190;10 126;11 70;12 28;13 12;14 5]; 
% % Somente os casos da subida da curva epidêmica 
%B=[1 3;2 8;3 28;4 75;5 221;6 291]; 
tempo=B(:,1);
tempo=tempo';
dados_obs=B(:,2);
dados_obs=dados_obs';

% Chute inicial para os parâmetros (ajuste conforme necessário)
chute_inicial = [1e-5, 1e-5];

% Usando a função fminsearch para minimizar a função de ajuste
params_otimizados = fminsearch(@(params) ajuste_SIR(params, tempo, dados_obs), chute_inicial);

% Parâmetros otimizados
beta_otimizado = params_otimizados(1);
gamma_otimizado = params_otimizados(2);

    S0 = 760;  % Exemplo de população inicial
    I0 = dados_obs(1);  % Valor inicial dos infectados
    R0 = 0;  % Nenhum recuperado inicialmente
   
    y0 = [S0; I0; R0];

% Resolvendo novamente as equações diferenciais com os parâmetros otimizados
 [t, y_ajustado] = ode45(@(t, y) modelo_SIR(t, y, beta_otimizado, gamma_otimizado), tempo, y0);

% =================== Visualizando os Resultados =========================
figure(1);
plot(tempo, dados_obs, 'o', t, y_ajustado(:,2), '-k');
legend('Dados Observados', 'Solução do Modelo SIR Ajustada');
xlabel('Tempo (dias)');
ylabel('Número de Infectados');

%% ===================== Funções Utilizadas ============================
% ======== Função de Ajuste para o Método dos Quadrados Mínimos ========
function sqe = ajuste_SIR(params, tempo, dados_obs)
    beta = params(1);
    gamma = params(2);
   
    % Condições iniciais
    S0 = 760;  % Exemplo de população inicial
    I0 = dados_obs(1);  % Valor inicial dos infectados
    R0 = 0;  % Nenhum recuperado inicialmente
   
    y0 = [S0; I0; R0];
   
    % Resolvendo as equações diferenciais usando ode45
    [~, y] = ode45(@(t, y) modelo_SIR(t, y, beta, gamma), tempo, y0);
   
    % Calculando a soma dos quadrados dos resíduos
    sqe = sum((y(:,2) - dados_obs').^2);
end

% =========== Função para Simular o Modelo SIR Sem Demografia ============
function dydt = modelo_SIR(~, y, beta, gamma)
    dydt = zeros(3,1);
    dydt(1) = -beta * y(1) * y(2);
    dydt(2) = beta * y(1) * y(2) - gamma * y(2);
    dydt(3) = gamma * y(2);
end
