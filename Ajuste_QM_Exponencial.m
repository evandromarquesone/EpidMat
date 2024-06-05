% Material Complementar do Livro: 
% Tópicos Matemáticos Aplicados na Modelagem em Epidemiologia - Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

%% =========== Exemplos de Aplicação do Método dos Quadrados Mínimos ======
% Ajuste dos casos de infectados de Influenza à função y=a*exp(b*x)
clc;  clear all;  close all; 

% % =========== Dados de Influenza no Internato ==========================
%B=[1 3;2 8;3 28;4 75;5 221;6 291;7 255;8 235;9 190;10 126;11 70;12 28;13 12;14 5]; 
%B=[1 3;2 8;3 28;4 75;5 221;6 291];
%B=[1 3;2 8;3 28;4 75;5 221]; 
%B=[1 3;2 8;3 28;4 75];
%B=[1 3;2 8;3 28];
B=[1 3;2 8];
xi=B(:,1);%vetor tempo. Primeira coluna da matriz B.
b=B(:,2);%vetor casos de indivíduos infectados. Segunda coluna da matriz B.

xi_dim=length(xi); %criando a variável com o tamanho do vetor xi
% ====================  Aplicando ln no vetor b ==========================
% Para alicar o método e determinar as constantes a e b de y=a*exp(b*x), é
% preciso transformar em uma equação linear. Aplicando ln em ambos os lados
% obtemos ln(y)=ln(a*exp(b*x))=ln(a)+b*x. Assim as funções são g1(x)=1 e 
% g2(x)=x. Para c=ln(a), temos uma equação
% linear: y=c+b*x. Precisamos então aplicar ln em cada elemento de b, e
% depois do processo, aplicar a exponencial, para finalmente ter os
% coeficientes.

lnb=zeros(xi_dim,1);
for i=1:xi_dim     
    lnb(i)=log(b(i));   
end

% =============== Gerando os Vetores ao Aplicar as =======================
g1=zeros(xi_dim,1);
g2=zeros(xi_dim,1);
for i=1:xi_dim
    g1(i)=1;
    g2(i)=xi(i);
end

% ========= Gerando a Matriz A com Produtos Internos das Funções==========
A=zeros(2);
A(1,1)=dot(g1,g1);
A(2,2)=dot(g2,g2);
A(1,2)=dot(g1,g2);
A(2,1)=dot(g2,g1);

% ======================== Gerando o vetor bBarra ========================
bBarra=zeros(2,1);
bBarra(1,1)=dot(lnb,g1); bBarra(2,1)=dot(lnb,g2); 
%% =========== Resolvendo o Sistema via Fatoração Cholesky ==============
%Para resolver o sistema com a fatoração Cholesky, baixe o programa 
% "solvespd" disponível em https://www.ime.unicamp.br/~pulino/ALESA/Matlab/ 
% E para entender melhor o uso das fatorações Cholesky e QR no método de Quadrados
% Mínimos, consulte: https://www.ime.unicamp.br/~marcia/AlgebraLinear/quadrados_minimos.html
 [G] = chol(A);
 [x] =solvespd(G,bBarra,2); 


% % ================== Construindo a Solução =============================
   xQMexpAlg=[exp(x(1)) x(2)]; 
% % ================== Plotando a Solução ==============================
figure(1)
x=0:0.1:xi_dim;
y= xQMexpAlg(1)*exp( xQMexpAlg(2)*x); %solução de y=ae^(bx) para dados ajustados
plot(x,y,'k')
hold on
plot(xi,b,'+ b')
axis([0 7 0 330])
% hold on
x2=xi_dim:0.1:6;
y2=xQMexpAlg(1)*exp( xQMexpAlg(2)*x2); %solução de y=ae^(bx) para além dos dados ajustados
plot(x2,y2,'-- k',3,28,'+ r',4,75,'+ r',5,221,'+ r',6, 291,'+ r')
%plot(x2,y2,'-- k')
%%B=[1 3;2 8;3 28;4 75;5 221;6 291];
legend('Solução de Quadrados Mínimos (QM)', 'Dados Observados Utilizados no Ajuste de QM',...
    'Solução de QM Estimando Casos Futuros','Dados Observados Não Utilizados no Ajuste')
xlabel('Tempo (dias)');
ylabel('Quantidade de Indivíduos Infectados');

%% ==== Plotando a Solução de Quadrados Mínimos com a Solução Anal[itica
% % == plotando junto a solução I(t)=Io*exp((beta*N-nu)*t)
figure(2)
x3=0:0.1:6;
n=763; %50 mil indivíduos
beta=1.66/763; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
nu=1/(2.2); %(a taxa é 1 dividido pelo tempo em que se fica infectado)
aux1=beta*n-nu;
y3=3*exp((aux1)*x3);
plot(x,y,'k',xi,b,'+ b',x3,y3,'r')
%plot(xi,b,'+ b')
axis([0 7 0 330])
hold on
axis([0 7 0 330])
xlabel('Tempo (dias)');
ylabel('Quantidade de Indivíduos Infectados');
legend('Solução de Quadrados Mínimos (QM)','Dados Observados Utilizados no Ajuste de QM','Solução Analítica')

%% ========= Construindo a Solução Ajustada e Calculando o Erro ==========
yi=zeros(xi_dim,1);
f1=zeros(xi_dim,1);
erro1=zeros(xi_dim,1);
erroQuad1=zeros(xi_dim,1);
for i=1:xi_dim
    yi(i)=xQMexpAlg(1)*( xQMexpAlg(2)).^xi(i); 
    erro1(i)=b(i)-yi(i);
    erroQuad1(i)=erro1(i)^2;
end
erroTotalF1=sum(erroQuad1);

% % ===================== Plotando os Erros ============================
msg = sprintf("Erro Quadrático Total Para Ajuste de %d Pontos", xi_dim);
 figure(3)
 plot(xi,erroQuad1)
xlabel('x'),
ylabel('Erro Quadrático da Solução de Quadrados Mínimos'),
title(msg),

% % ============= Mensagens Exibidas no Display =========================

disp('O erro total quadrático é:')
disp(erroTotalF1) 

disp('O valor aproximado de infectados para o terceiro dia é:')
disp(xQMexpAlg(1)*exp( xQMexpAlg(2)*3))

disp('O valor aproximado de infectados para o quarto dia é:')
disp(xQMexpAlg(1)*exp( xQMexpAlg(2)*4))

disp('O valor aproximado de infectados para o quinto dia é:')
disp(xQMexpAlg(1)*exp( xQMexpAlg(2)*5))

disp('O valor aproximado de infectados para o sexto dia é:')
disp(xQMexpAlg(1)*exp( xQMexpAlg(2)*6))

disp('O valor real de infectados no sexto dia é 291.')
