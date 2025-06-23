%% ==== Programa de Simulação Estocástica do Modelo Sir com Demografia ====
% o modelo está disponível no arquivo sir_ComDemografia.m
clc;  clear all;  close all;

% =================== Condições Iniciais do Programa ====================
NInt=40; %número de iterações
T=1000000;%tempo em dias
Int=1000000; %intervalo de tempo, recomendável ser igual a T 
dt=T/Int;
cont=0; 
n=5000; %5 mil indivíduos
%% =========================== Parâmetros ================================

% Lembrando os valores e significados de cada parâmetro
% n=5000; %50 mil indivíduos 
% beta=0.1/n; %taxa per-capita de 10% dos encontros tornarem-se contaminados)
% mu=1/(80*365); %(expectativa de 80 anos de vida)
% nu=1/60; %(a taxa é 1 dividido pelo tempo em que se fica infectado)


% ========== vetor com os valores dos parâmetros =========================
% %            muN,      beta       mu        nu           
vTaxaPadrao=[n/(80*365); 0.1/n; 1/(80*365); 1/60];
muN=vTaxaPadrao(1,1); beta=vTaxaPadrao(2,1); mu=vTaxaPadrao(3,1); 
nu=vTaxaPadrao(4,1); 

% ============  Declarando matrizes nulas do problema ====================
Tstop=zeros(NInt,1);
S=zeros(T+1,NInt);
I=zeros(T+1,NInt);
R=zeros(T+1,NInt);


C=zeros(T+1,NInt);
Sn1=zeros(T+1,1); Sn2=zeros(T+1,1); Sn3=zeros(T+1,1);
In1=zeros(T+1,1); In2=zeros(T+1,1); In3=zeros(T+1,1);
Rn1=zeros(T+1,1); Rn2=zeros(T+1,1); Rn3=zeros(T+1,1);
%Yn1=zeros(T+1,NInt); Yn2=zeros(T+1,NInt); Yn3=zeros(T+1,NInt);
Cn1=zeros(T+1,NInt); Cn2=zeros(T+1,NInt); Cn3=zeros(T+1,NInt);

SM=zeros(NInt+1,1);IM=zeros(NInt+1,1); RM=zeros(NInt+1,1); 

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
     % %==== condição inicial Equilíbrio Endêmico
%      L(1,i)=fix((fi/(g2+muC))*((Rg-1)/ErreZeroD));
%      S(1,i)=(fi-(g2+muC)*L(1,i))/(muC+nu); 
%      Ad(1,i)=(g2/(g3+g4+muC))*L(1,i) ;
%      I(1,i)=(g2*(g3*(g6+theta2+muC)+g4*g6)*L(1,i))/((g3+g4+muC)*((g5+theta1+delta+muC)*(g6+theta2+muC)-g5*g6));
%      Ac(1,i)=((g2*(g4*(g5+theta1+delta+muC)+g3*g5))/((g3+g4+muC)*((g5+theta1+delta+muC)*(g6+theta2+muC)-g5*g6)))*L(1,i);
     %Y=(chiI*n*L*(g2+muC))/((fi/(muC+nu))*alphaM*n+chiI*(g2+muC)*L);
     
     % condição inicial S(0)=4999, I(0)=1, R(0)=0 
     S(1,i)=4999; 
     I(1,i)=1;
     R(1,i)=0;


 %================== Armazenando as Informações ==========================
    TT=0; %tempo de mudança/transição
    posant=1; %posição anterior
    Sant=S(1,i);
    Iant=I(1,i);
    Rant=R(1,i);
    Cant=Sant+Iant+Rant; %população de indivíduos

%============================== Taxas ====================================
disp(i)

    while TT<=T
    r1=muN; r2=mu*Sant; %entrada e saída natural de S
    r3=beta*Iant*Sant; %saída de S por contágio 
    r4=nu*Iant;r5=mu*Iant; %saídas de I
    r6=mu*Rant; %saída de R
                   
    Resultante=r1+r2+r3+r4+r5+r6;
    %===== Cálculo do Tempo de Transição/Infectados + Assintomáticos = 0 =====
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
     
    %======================== Probabilidades =================================
        Pr1=r1/Resultante;   Pr2=r2/Resultante;  Pr3=r3/Resultante;
        Pr4=r4/Resultante;   Pr5=r5/Resultante;  Pr6=r6/Resultante;

        tr=rand(1);

        if tr<Pr1
            Sant=Sant+1;
        elseif tr<Pr1+Pr2
            Sant=Sant-1;
        elseif tr<Pr1+Pr2+Pr3
            Sant=Sant-1;
            Iant=Iant+1;
        elseif tr<Pr1+Pr2+Pr3+Pr4
            Iant=Iant-1; 
        elseif tr<Pr1+Pr2+Pr3+Pr4+Pr5
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
    
end

% Comandos para salvar os dados de simulação e todas populações em todas
% iterações.

% txtfile=['parametro' '.txt'];
% f=fopen(txtfile, 'w');
% fprintf(f,'Pop.\t ni\t beta\t gama\t sigma\t I\t T\t cont\n');
% fclose(f);
% parametros=[Np niest b g s I T cont];
% save (txtfile, 'parametros', '-ascii','-append');
% 
% txtfile=['Pops' '.txt'];
% g=fopen(txtfile, 'w');
% fprintf(g,'X\t Hnv\t Hv\t Ynv\t Yv\t  Znv\t Zv\t V\t Va\t D\t Da\n');
% % Hnv\t Hv\t Ynv\t Yv\t  Znv\t Zv\t V\t Va\t D\t Da\n'
% fclose(g);
% pops=[X Hnv Hv Ynv Yv Znv Zv V Va D Da];
% save (txtfile, 'pops', '-ascii','-append');
%
%=========  Montagem do Vetor Tempo ===================
for k=1:Int+1%1001
    temp(k)=(k-1)*dt;
end
% 
% soma=HM+YM; %soma de exposto com infectado
% Tstop1=temp(find(soma<1,1));
% Tstop05=temp(find(soma<0.5,1));
% Tstop01=temp(find(soma<0.1,1));
% %plot(temp,XM,temp,XMed,temp,Xq,temp,XQ)

% ======================== Plotagem de Gráficos ========================= 
figure(1)
title('Media, mediana, quantis 0.025 e 0.975 com tempo continuo')
subplot(2,2,1), plot(temp,SM,temp,SMed,temp,Sq,temp,SQ)
title('Media, mediana, quantis 0.025 e 0.975 com tempo continuo')
xlabel('tempo (dias)')
ylabel('Suscetiveis')
subplot(2,2,2), plot(temp,IM,temp,IMed,temp,Iq,temp,IQ)
xlabel('tempo (dias)')
ylabel('Infectados')
subplot(2,2,3), plot(temp,RM,temp,RMed,temp,Rq,temp,RQ)
xlabel('tempo (dias)')
ylabel('Recuperados')
subplot(2,2,4), plot(temp,CM,temp,CMed,temp,Cq,temp,CQ)
xlabel('tempo (dias)')
ylabel('tamanho pop.')


figure(2)
title('Trajetorias aleatorias tempo continuo')
subplot(2,2,1), plot(temp,Sn1,temp,Sn2,temp,Sn3)
title('Trajetorias aleatorias tempo continuo' )
xlabel('tempo (dias)')
ylabel('Suscetiveis')
subplot(2,2,2), plot(temp,In1,temp,In2,temp,In3)
xlabel('tempo (dias)')
ylabel('Infectados')
subplot(2,2,3), plot(temp,Rn1,temp,Rn2,temp,Rn3)
xlabel('tempo (dias)')
ylabel('Recuperados')
subplot(2,2,4), plot(temp,Cn1,temp,Cn2,temp,Cn3)
xlabel('tempo (dias)')
ylabel('tamanho pop.')
% subplot(3,3,9), plot(temp,Nn1,temp,Nn2,temp,Nn3)
% xlabel('tempo (dias)')
% ylabel('tamanho pop.')
% subplot(3,3,5), plot(temp,Dn1,temp,Dn2,temp,Dn3)
% xlabel('tempo (dias)')
% ylabel('Desperdicio')
% subplot(3,3,7), plot(temp,Van1,temp,Van2,temp,Van3)
% xlabel('tempo (dias)')
% ylabel('Vacinados acumulados ex-susc.')
% subplot(3,3,8), plot(temp,Dan1,temp,Dan2,temp,Dan3)
% xlabel('tempo (dias)')
% ylabel('Desperdicio acumulados')
