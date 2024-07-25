%% Material Complementar do Livro: 
% Introdução à Epidemiologia Matemática: Métodos em Estudos Transversais

% = Outros Materiais estão disponíveis em https://linktr.ee/livroepidmat =

%% ==== Programa de Simulação Estocástica do Modelo SIR com Demografia ====
% Código implementado de acordo com a descrição do algoritmo apresentada na
% página 88 do livro.

 clc;  clear all;  close all;

% =================== Condições Iniciais do Programa ====================
NInt=40; %número de iterações
T=100;%tempo em dias
Int=T; %intervalo de tempo, recomendável ser igual a T 
dt=T/Int;
cont=0; 
CI=50000; %população total inicial de indivíduos
% =========================== Parâmetros ================================
% Com fi=50 000/(80*365), mu=1/(80*365), beta=0,1/50 000, nu=1/60 tem-se Ro=5,9877 e chi=0,167; 
                                                                                       
%            fi=mu*N,         mu,     beta,   nu,     
vTaxaPadrao=[CI/(80*365); 1/(80*365);0.1/CI; 1/60];
fi=vTaxaPadrao(1,1); mu=vTaxaPadrao(2,1); beta=vTaxaPadrao(3,1); nu=vTaxaPadrao(4,1);

% ====================== Definindo Ro e Chi===========================
Ro=(beta*CI)/(nu+mu);
chi=1/Ro;

% ============  Declarando matrizes nulas do problema ====================
Tstop=zeros(NInt,1);
S=zeros(T+1,NInt);
I=zeros(T+1,NInt);
R=zeros(T+1,NInt);


C=zeros(T+1,NInt);
Sn1=zeros(T+1,1); Sn2=zeros(T+1,1); Sn3=zeros(T+1,1);
In1=zeros(T+1,1); In2=zeros(T+1,1); In3=zeros(T+1,1);
Rn1=zeros(T+1,1); Rn2=zeros(T+1,1); Rn3=zeros(T+1,1);
Cn1=zeros(T+1,NInt); Cn2=zeros(T+1,NInt); Cn3=zeros(T+1,NInt);

SM=zeros(NInt+1,1);IM=zeros(NInt+1,1);RM=zeros(NInt+1,1); 

Sdp=zeros(NInt+1,1);Idp=zeros(NInt+1,1);Rdp=zeros(NInt+1,1);

Sq=zeros(NInt+1,1);Iq=zeros(NInt+1,1);Rq=zeros(NInt+1,1);

SQ=zeros(NInt+1,1);IQ=zeros(NInt+1,1);RQ=zeros(NInt+1,1);

SMed=zeros(NInt+1,1);IMed=zeros(NInt+1,1);RMed=zeros(NInt+1,1);

CM=zeros(NInt+1,1);Cdp=zeros(NInt+1,1);Cq=zeros(NInt+1,1);CQ=zeros(NInt+1,1);
CMed=zeros(NInt+1,1);

temp=zeros(NInt+1,1);

% ============== Gerando Números Aleatórios ==============================
n1=randi(NInt);
n2=randi(NInt-1);
n3=randi(NInt-2);
    if n2>=n1
        n2=n2+1;
    end
    if n3>=min(n1,n2)
        n3=n3+1;
    elseif n3>=max(n1,n2)
        n3=n3+1;
    end
    
 %= Rodando o Programa a Partir da Soluçao Endêmica como Condição Inicial =
 for i=1:NInt     
     % condição inicial S(0)=49 999, I(0)=1 e R(0)=0 
     S(1,i)=49999;    
     I(1,i)=1;
     R(1,i)=0;
 %================== Armazenando as Informações ==========================
    TT=0; %tempo de mudança/transição
    posant=1; %posição anterior
    Sant=S(1,i);
    Iant=I(1,i);
    Rant=R(1,i);

    Cant=Sant+Iant+Rant; %população total

%===================== Imprimindo no Display =============================
disp(i)

% %PopInicialNames={'S';'L';'Ad';'I';'R'};
% CondInicialName={'ValorInicial'};
% %S=Sant; L=Lant; Ad=Adant;I=Iant;R=Rant;Total=Cant;
% Si=S(1,i); Li=L(1,i); Adi=Ad(1,i);Ii=I(1,i);Ri=R(1,i);Totali=Cant;
% table(Si,Li,Adi,Ii,Ri,Totali,'RowNames',CondInicialName)

%============================== Taxas ====================================
    while TT<=T
    %r0=mu*Cant;
    %r1=fi %original    e r1=fi*Cant modificado
    r1=fi; r2=mu*Sant; %entradas e saída naturais de S
    r3=beta*Sant*Iant; %saída de S por contágio
    r4=mu*Iant;r5=nu*Iant; %saídas de I
    r6=mu*Rant; %saída de R

% somando os eventos para depois calcular as propabilidades                   
    Resultante=r1+r2+r3+r4+r5+r6;
    vetorTaxas=[r1;r2;r3;r4;r5;r6;Resultante];

%     % ===== Cálculo do Tempo de Transição/Infectados = 0 =====
     TM=exprnd(1/Resultante);%tempo de mudança/transição
         if Iant+Rant>0 % do problema
             Tstop(i,1)=TT;
         end
         TT=TT+TM; %correspondente ao algoritmo de Gillespie 
         pos=ceil(TT/dt)+1;
         if pos>T
             pos=Int+1;
         end
         Sant;

         S(posant:pos,i)=Sant;        
         I(posant:pos,i)=Iant;
         R(posant:pos,i)=Rant;
     
    %======================== Probabilidades ==============================
        Pr1=r1/Resultante;   Pr2=r2/Resultante;  Pr3=r3/Resultante;
        Pr4=r4/Resultante;   Pr5=r5/Resultante;  Pr6=r6/Resultante;

       probTotal=Pr1+Pr2+Pr3+Pr4+Pr5+Pr6;
      vetorProb=[Pr1;Pr2;Pr3;Pr4;Pr5;Pr6;probTotal];

    % =============== Iniciando o Sistema Dinâmico =====================        
        tr=rand(1);

        if tr<Pr1
            Sant=Sant+1;
        elseif tr<Pr1+Pr2
            Sant=Sant-1;
        elseif tr<Pr1+Pr2+Pr3
            Sant=Sant-1;
            Iant=Iant+1;
        elseif tr<(Pr1+Pr2+Pr3+Pr4)
            Iant=Iant-1;
        elseif tr<(Pr1+Pr2+Pr3+Pr4+Pr5)
            Iant=Iant-1;
            Rant=Rant+1;           
        else 
            Rant=Rant-1;
         end    
        Cant=Sant+Iant+Rant;
        posant=pos;
    end
     

%========== Contando Quantas vezes Tstop ficou menor do que T ============
     if Tstop(i,1)<T
         cont=cont+1;
     end   

 end
 
 % =========== Imprimindo Equilíbrio Determinístico e Ro =================
     Seq=(CI/Ro);
     Ieq=(mu/beta)*(Ro-1);
     Req=(nu/beta)*(Ro-1);
     TotalEq=Seq+Ieq+Req;
   CondFinalName={'EqDet'};
   table(Seq,Ieq,Req,TotalEq,'RowNames',CondFinalName)
   
   RoNames={'ValorAbsoluto'};
  table(Ro,chi,'RowNames',RoNames)

 % ==== Organizando para formar as soluções aleatórias selecionadas ======
 
Sn1=S(:,n1);           Sn2=S(:,n2);          Sn3=S(:,n3);
In1=I(:,n1);           In2=I(:,n2);          In3=I(:,n3);
Rn1=R(:,n1);           Rn2=R(:,n2);          Rn3=R(:,n3);

Cn1=Sn1+In1+Rn1;  Cn2=Sn2+In2+Rn2;  Cn3=Sn3+In3+Rn3;

C=S+I+R;

%=========== Cálculo das Médias, Medianas, Desvios Padrões ============
for k=1:Int+1
    %Médias
    SM(k)=mean(S(k,:));
    IM(k)=mean(I(k,:));
    RM(k)=mean(R(k,:));

    CM(k)=mean(C(k,:));
    
    %Desvio Padrão
    Sdp(k)=std(S(k,:));
    Idp(k)=std(I(k,:));
    Rdp(k)=std(R(k,:));

    Cdp(k)=std(C(k,:));
    
    %Quantis
    sor=sort(S(k,:)); %ordena a linha em ordem crescente
    Sq(k)=(sor(0.025*NInt)+sor(0.025*NInt+1))/2;
    SQ(k)=(sor(0.975*NInt)+sor(0.975*NInt+1))/2;
    SMed(k)=(sor(0.5*NInt)+sor(0.5*NInt+1))/2;    
      
    sor=sort(I(k,:)); 
    Iq(k)=(sor(0.025*NInt)+sor(0.025*NInt+1))/2;
    IQ(k)=(sor(0.975*NInt)+sor(0.975*NInt+1))/2;
    IMed(k)=(sor(0.5*NInt)+sor(0.5*NInt+1))/2;
    
    sor=sort(R(k,:)); 
    Rq(k)=(sor(0.025*NInt)+sor(0.025*NInt+1))/2;
    RQ(k)=(sor(0.975*NInt)+sor(0.975*NInt+1))/2;
    RMed(k)=(sor(0.5*NInt)+sor(0.5*NInt+1))/2;
    
    
    sor=sort(C(k,:)); 
    Cq(k)=(sor(0.025*NInt)+sor(0.025*NInt+1))/2;
    CQ(k)=(sor(0.975*NInt)+sor(0.975*NInt+1))/2;
    CMed(k)=(sor(0.5*NInt)+sor(0.5*NInt+1))/2;
    
    %=========  Montagem do Vetor Tempo ===================
%for k=1:Int+1%1001
    temp(k)=(k-1)*dt;
%end
end

%% ======================== Plotagem de Gráficos ========================= 
figure(1)
title('Media, mediana, quantis 0.025 e 0.975 com tempo continuo')
subplot(3,1,1), plot(temp,SM,temp,SMed,temp,Sq,temp,SQ)
title('Media, mediana, quantis 0.025 e 0.975 com tempo continuo')
xlabel('tempo (dias)')
ylabel('Suscetiveis')
subplot(3,1,2), plot(temp,IM,temp,IMed,temp,Iq,temp,IQ)
xlabel('tempo (dias)')
ylabel('Infectados')
subplot(3,1,3), plot(temp,RM,temp,RMed,temp,Rq,temp,RQ)
xlabel('tempo (dias)')
ylabel('Recuperados')

figure(2)
title('Trajetorias aleatorias no tempo continuo')
subplot(3,1,1), plot(temp,Sn1,temp,Sn2,temp,Sn3)
title('Trajetorias aleatorias tempo continuo' )
xlabel('tempo (dias)')
ylabel('Suscetiveis')
subplot(3,1,2), plot(temp,In1,temp,In2,temp,In3)
xlabel('tempo (dias)')
ylabel('Infectados')
subplot(3,1,3), plot(temp,Rn1,temp,Rn2,temp,Rn3)
xlabel('tempo (dias)')
ylabel('Recuperados')

%% =========== Organizando e Salvando os dados ==========================
saidasSsemCon=[SM(Int+1,1);SMed(Int+1,1);Sq(Int+1,1);SQ(Int+1,1);Ro];
saidasCsemCon=[CM(Int+1,1);CMed(Int+1,1);Cq(Int+1,1);CQ(Int+1,1)];
 

 % txtfile=['saidasSsemCon' '.txt'];
 % save(txtfile, 'saidasSsemCon', '-ascii','-append');
 % 
 % txtfile=['saidasCsemCon' '.txt'];
 % save(txtfile, 'saidasCsemCon', '-ascii','-append');



Matrizq02=[Sq Iq Rq Cq];
MatrizQ97=[SQ IQ RQ CQ];
MatrizM=[SM IM RM CM];
MatrizMed=[SMed IMed RMed CMed];

%% ======== Comandos para Salvar as Matrizes em Arquivos .txt ============
% txtfile=['Matrizq02' '.txt'];
%  save(txtfile, 'Matrizq02', '-ascii','-append');
% 
%  txtfile=['MatrizQ97' '.txt'];
%  save(txtfile, 'MatrizQ97', '-ascii','-append');
% 
%  txtfile=['MatrizM' '.txt'];
%  save(txtfile, 'MatrizM', '-ascii','-append');
% 
%  txtfile=['MatrizMed' '.txt'];
%  save(txtfile, 'MatrizMed', '-ascii','-append');
