% Material Complementar do Livro: 
% Tópicos Matemáticos Aplicados na Modelagem em Epidemiologia - Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

%% =========== Exemplos de Aplicação do Método dos Quadrados Mínimos ======
% % ======= Ajuste de Quadrados Mínimos para a função y=ax^2+bx+c =========
clc;  clear all;  close all; 

% % =========== Dados de Influenza no Internato ==========================
%B=[1 3;2 8;3 28;4 75;5 221;6 291;7 255;8 235;9 190;10 126;11 70;12 28;13 12;14 5]; 
%B=[1 3;2 8;3 28;4 75;5 221;6 291];
%B=[1 3;2 8;3 28;4 75;5 221]; 
%B=[1 3;2 8;3 28;4 75];
B=[1 3;2 8;3 28];
%B=[1 3;2 8]; %com esse vetor, a Matriz A é singular
xi=B(:,1);%vetor tempo. Primeira coluna da matriz B.
b=B(:,2);%vetor casos de indivíduos infectados. Segunda coluna da matriz B.

xi_dim=length(xi); %criando a variável com o tamanho do vetor xi
% ====================  Aplicando ln no vetor b ==========================
% Para alicar o método e determinar as constantes a e b de y=x^a, é
% preciso transformar em uma equação linear. Aplicando ln em ambos os lados
% obtemos ln(y)=ln(x^a)=a*ln(x). Assim a única função é g1(x)=1. 
% Precisamos então aplicar ln em cada elemento de b e de x, e
% depois do processo, aplicar a exponencial, para finalmente ter a solução.
% lnb=zeros(xi_dim,1);
% for i=1:xi_dim     
%     lnb(i)=log(b(i)); 
% end

% =============== Gerando os Vetores ao Aplicar as =======================
g1=zeros(xi_dim,1);
g2=zeros(xi_dim,1);
g3=zeros(xi_dim,1);
for i=1:xi_dim
    g1(i)=xi(i)^2;
    g2(i)=xi(i);
    g3(i)=1;
end

% ========= Gerando a Matriz A com Produtos Internos das Funções==========
A=zeros(3);
A(1,1)=dot(g1,g1); A(2,2)=dot(g2,g2); A(3,3)=dot(g3,g3);
A(1,2)=dot(g1,g2); A(2,1)=A(1,2);
A(1,3)=dot(g1,g3); A(3,1)=A(1,3);
A(2,3)=dot(g2,g3); A(3,2)=A(2,3);

% ======================== Gerando o vetor bBarra ========================
bBarra=zeros(3,1);
bBarra(1,1)=dot(b,g1); bBarra(2,1)=dot(b,g2); bBarra(3,1)=dot(b,g3); 

%% ===================== Resolvendo o Sistema Linear =====================
% ===================== Solução via Fatoração Cholesky ===================
% % ========= A matriz A do sistema deve ser positiva definida ===========
% %Para resolver o sistema com a fatoração Cholesky, baixe o programa 
% % "solvespd" disponível em https://www.ime.unicamp.br/~pulino/ALESA/Matlab/ 
 [G] = chol(A);
 [x_Sol] =solvespd(G,bBarra,3); 
% %=================== Solução via Fatoração QR ===========================
% %Para resolver o sistema com a fatoração QR, baixe o programa 
% % "trisuplin" disponível em https://www.ime.unicamp.br/~pulino/ALESA/Matlab/
% [Q,R]=qr(A); %determinando a fatoração QR de A
% bBarra_aux=zeros(3,1);
% bBarra_aux=Q'*bBarra;
% [x_Sol] = trisuplin(R,bBarra_aux);
% [x_Sol]=A\bBarra;

% % ================== Construindo a Solução =============================
   xQM_Coef=[x_Sol(1) x_Sol(2) x_Sol(3)]; 
% % ================== Plotando a Solução ==============================
figure(1)
x=0:0.1:xi_dim;
y= xQM_Coef(1)*x.*x+xQM_Coef(2)*x+xQM_Coef(3); %solução de y=ax^2+bx+c para dados ajustados
plot(x,y,'k')
hold on
plot(xi,b,'+ b')
axis([0 7 0 330])
% hold on
x2=xi_dim:0.1:7;
y2= xQM_Coef(1)*x2.*x2+xQM_Coef(2)*x2+xQM_Coef(3); %solução de y=ax^2+bx+c para dados ajustados
%plot(x2,y2,'-- k',3,28,'+ r',4,75,'+ r',5,221,'+ r',6, 291,'+ r')
plot(x2,y2,'-- k',4,75,'+ r',5,221,'+ r',6, 291,'+ r')
% plot(x2,y2,'-- k',5,221,'+ r',6, 291,'+ r')
%plot(x2,y2,'-- k',6, 291,'+ r')
plot(x2,y2,'-- k')
%%B=[1 3;2 8;3 28;4 75;5 221;6 291];
legend('Solução de Quadrados Mínimos (QM)', 'Dados Observados Utilizados no Ajuste de QM',...
    'Solução de QM Estimando Casos Futuros','Dados Observados Não Utilizados no Ajuste')
xlabel('Tempo (dias)');
ylabel('Quantidade de Indivíduos Infectados');

%% ==== Plotando a Solução de Quadrados Mínimos com a Solução Analítica
% % == plotando junto a solução I(t)=Io*exp((beta*N-nu)*t)
figure(2)
x3=0:0.1:7;
n=763; %50 mil indivíduos
beta=1.66/763; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
nu=1/(2.2); %(a taxa é 1 dividido pelo tempo em que se fica infectado)
aux1=beta*n-nu;
y3=3*exp((aux1)*x3);
plot(x,y,'k',xi,b,'+ b',x3,y3,'r')
%plot(xi,b,'+ b')
axis([0 7 0 330])
hold on  % para plotar a solução de QM polinomial com a exponencial para
          % 6 pontos
y4=1.3078*exp(0.9661*x3);
plot(x3,y4,'- b')
axis([0 7 0 330])
% hold on %para obter a solução média, descomente
% y5= xQM_Coef(1)*x3.*x3+xQM_Coef(2)*x3+xQM_Coef(3); %solução de y=ax^2+bx+c para dados ajustados
% y6=(y4+y5)/2;
% plot(x3,y6,'- g')
xlabel('Tempo (dias)');
ylabel('Quantidade de Indivíduos Infectados');
%legend('Solução de Quadrados Mínimos (QM)','Dados Observados Utilizados no Ajuste de QM','Solução Analítica')
 legend('Solução de Quadrados Mínimos (QM) Polinomial','Dados Observados Utilizados no Ajuste de QM','Solução Analítica', 'Solução de QM Exponencial')
%legend('Solução de Quadrados Mínimos (QM) Polinomial','Dados Observados Utilizados no Ajuste de QM','Solução Analítica', 'Solução de QM Exponencial','Solução Média de QM')

%% ========= Construindo a Solução Ajustada e Calculando o Erro ==========
yi=zeros(xi_dim,1);
f1=zeros(xi_dim,1);
erro1=zeros(xi_dim,1);
erroQuad1=zeros(xi_dim,1);
for i=1:xi_dim
   % yi(i)=xQMexpAlg(1)*( xQMexpAlg(2)).^xi(i); 
   y(i)= xQM_Coef(1)*x(i)^2+xQM_Coef(2)*x(i)+xQM_Coef(3); %solução de y=ax^2+bx+c para dados ajustados
    erro1(i)=b(i)-yi(i);
    erroQuad1(i)=erro1(i)^2;
end
erroTotalF1=sum(erroQuad1);

% % ===================== Plotando os Erros ============================
msg = sprintf("Erro Quadrático Total Para Ajuste de %d Pontos", xi_dim);
% mensagem a ser exibida na legenda, de acordo com a quantidade de pontos
% utilizada no ajuste
figure(3)
 plot(xi,erroQuad1)
xlabel('x'),
ylabel('Erro Quadrático da Solução de Quadrados Mínimos'),
title(msg),


% % ============= Mensagens Exibidas no Display =========================

disp('O erro total quadrático é:')
disp(erroTotalF1) 

disp('O valor aproximado de infectados para o terceiro dia é:')
disp(xQM_Coef(1)*3^2+xQM_Coef(2)*3+xQM_Coef(3))

disp('O valor aproximado de infectados para o quarto dia é:')
disp(xQM_Coef(1)*4^2+xQM_Coef(2)*4+xQM_Coef(3))

disp('O valor aproximado de infectados para o quinto dia é:')
disp(xQM_Coef(1)*5^2+xQM_Coef(2)*5+xQM_Coef(3))

disp('O valor aproximado de infectados para o sexto dia é:')
disp(xQM_Coef(1)*6^2+xQM_Coef(2)*6+xQM_Coef(3))

%disp('O valor real de infectados no sexto dia é 291.')
