%{

    Code for reproduced paper results: COVID-19 Analisis With The Generalized SEIR Model Implemented With Rungge Kutta of 4th order (RK4)    
    S --> Susceptible
    P --> Insusceptible
    E --> Exposed
    I --> Infected
    Q --> Quarantined
    R --> Recovered
    D --> Death
    
    References:
    https://arxiv.org/pdf/2002.06563.pdf
%}


function generalizedSeirModel()
   

    close all
    clear all
    tMax = 55;    
    h = tMax/10000; %step size
    
    %Initial conditions............
    N = 8000000;
    E = 2;
    I = 10;
    S = N - E - I;
    Q = 0;
    R = 0;
    D = 0;
    P = 0;     
    % Some data is stimated in [1]
    alpha = 0.085; %Protection Rate
    betha = 1; %Infection Rate
    gamma = 1/2; % gamma^-1 Latent time
    delta = 1/(7.4) ;% delta^-1 Quarantine time in Wuhan
    kappa = 0.01; % Mortality rate
    lamda = 0.01; %Curate rate
    count = 1;
    
    
% slo ===> Matrix (j x 7) that contained the slopes   
%     k ---> (j,1)
%     l ---> (j,2)
%     n ---> (j,3)
%     m ---> (j,4)
%     w ---> (j,5)
%     z ---> (j,6)
%     r ---> (j,7)
% where, j  = 1,2,3,4
    
    for i = 1 : h : tMax-1
        
        vS(count) = S;
        vE(count) = E;
        vI(count) = I;
        vQ(count) = Q;
        vR(count) = R;
        vD(count) = D;
        vP(count) = P;
       
        slo = slopesComputation(betha,alpha,gamma,delta,lamda,kappa,S,E,I,Q,N,h);
        
        Snew = S + 1/6 *(slo(1,1) + 2*slo(2,1) + 2*slo(3,1) + slo(4,1)).*h;
        Enew = E + 1/6 *(slo(1,2) + 2*slo(2,2)+ 2*slo(3,2) + slo(4,2)).*h;
        Inew = I + 1/6 *(slo(1,3) + 2*slo(2,3) + 2*slo(3,3) + slo(4,3)).*h;
        Qnew = Q + 1/6 *(slo(1,4) + 2*slo(2,4) + 2*slo(3,4) + slo(4,4)).*h;
        Rnew = R + 1/6 *(slo(1,5) + 2*slo(2,5) + 2*slo(3,5) + slo(4,5)).*h;
        Dnew = D + 1/6 *(slo(1,6) + 2*slo(2,6) + 2*slo(3,6) + slo(4,6)).*h;
        Pnew = P + 1/6 *(slo(1,7) + 2*slo(2,7) + 2*slo(3,7) + slo(4,7)).*h;
        
        S = Snew;
        E = Enew;
        I = Inew;
        Q = Qnew;
        R = Rnew;
        D = Dnew;
        P = Pnew;
        
        
        count = count + 1;
        
    end
  
    
    
    ploting(vS,vE,vI,vQ,vR,vD,vP,tMax);

end


function slope = slopesComputation(betha,alpha,gamma,delta,lamda,kappa,S,E,I,Q,N,h)
    
    k = zeros(1,4);
    l = zeros(1,4);
    n = zeros(1,4);
    m = zeros(1,4);
    w = zeros(1,4);
    z = zeros(1,4);
    r = zeros(1,4);
    slopes = zeros(4,7);
    
    k(1) = (-1*betha/N)*S*I - alpha*S;
    l(1) = (betha/N)*S*I - gamma*E;
    n(1) = gamma*E - delta*I;
    m(1) = delta*I - lamda*Q - kappa*Q;
    w(1) = lamda*Q;
    z(1) = kappa*Q;
    r(1) = alpha*S;
    
    k(2) = (-1*betha/N) * ( S + k(1)*(h/2) ) * (I + n(1) *(h/2)) - alpha*(S + k(1)*(h/2) );
    l(2) =  (betha/N) * ( S + k(1)*(h/2) ) * (I + n(1) *(h/2)) - gamma * (E + l(1)*(h/2) );
    n(2) = gamma * (E + l(1)*(h/2) ) - delta * (I + n(1) * (h/2));
    m(2) = delta * (I + n(1) * (h/2)) - lamda * (Q + m(1) *(h/2)) - kappa*(Q + m(1) *(h/2));
    w(2) = lamda*(Q + m(1) *(h/2));
    z(2) = kappa*(Q + m(1) *(h/2));
    r(2) = alpha*( S + k(1)*(h/2) );
    
    k(3) = (-1*betha/N) * ( S + k(2)*(h/2) ) * (I + n(2) *(h/2)) - alpha*(S + k(2)*(h/2) );
    l(3) =  (betha/N) * ( S + k(2)*(h/2) ) * (I + n(2) *(h/2)) - gamma * (E + l(2)*(h/2) );
    n(3) = gamma * (E + l(2)*(h/2) ) - delta * (I + n(2) * (h/2));
    m(3) = delta * (I + n(2) * (h/2)) - lamda * (Q + m(2) *(h/2)) - kappa*(Q + m(2) *(h/2));
    w(3) = lamda*(Q + m(2) *(h/2));
    z(3) = kappa*(Q + m(2) *(h/2));
    r(3) = alpha*( S + k(2)*(h/2) );
    
    k(4) = (-1*betha/N) * ( S + k(3)*(h) ) * (I + n(3)*(h)) - alpha*(S + k(3)*(h) );
    l(4) =  (betha/N) * ( S + k(3)*(h) ) * (I + n(3) *(h)) - gamma * (E + l(3)*(h) );
    n(4) = gamma * (E + l(3)*(h) ) - delta * (I + n(3) * (h));
    m(4) = delta * (I + n(3) * (h)) - lamda * (Q + m(3) *(h)) - kappa*(Q + m(3) *(h));
    w(4) = lamda*(Q + m(3) *(h));
    z(4) = kappa*(Q + m(3) *(h));
    r(4) = alpha*( S + k(3)*(h) );
    
    slopes(:,1) = k;
    slopes(:,2) = l;
    slopes(:,3) = n;
    slopes(:,4) = m;
    slopes(:,5) = w;
    slopes(:,6) = z;
    slopes(:,7) = r;
    
    
    slope = slopes;  

end


function ploting(vS,vE,vI,vQ,vR,vD,vP,tMax)

    t = linspace(0,tMax,length(vS));
    matrix = zeros(7,length(vS));
    matrix(1,:) = vS;
    matrix(2,:) = vE;
    matrix(3,:) = vI;
    matrix(4,:) = vQ;
    matrix(5,:) = vR;
    matrix(6,:) = vD;
    matrix(7,:) = vP;
    
    for i = 1 : 7
        
        figure(i)
        plot(t,matrix(i,:),'b')
        xlabel('$\textbf{Time(days)}$','interpreter','latex')
        ylabel('$\textbf{Number of people}$','interpreter','latex')
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin'; 
        
        
        
       switch i         
           case 1           
              title("PERSONAS SUCEPTIBLES COVID-19");  
           case 2
               title("PERSONAS EXPUESTAS COVID-19"); 
           case 3
               title("PERSONAS INFECTADAS COVID-19");
           case 4
               title("PERSONAS EN CUARENTENA COVID-19");
           case 5
               title("PERSONAS RECUPERADAS COVID-19");
           case 6 
               title("PERSONAS MUERTAS COVID-19")
           case 7
               title("PERSONAS AISLADAS COVID-19");
               
        end
        hold on 
        
    end
    
    figure(8)
    plot(t,matrix(2,:) + matrix(3,:) ,'r')
    title("PERSONAS EXPUESTAS + INFECTADAS COVID-19");
    xlabel('$\textbf{Time(days)}$','interpreter','latex')
    ylabel('$\textbf{Number of people}$','interpreter','latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin'; 
    hold on
    plot(t,matrix(4,:),'b');
    legend("Expuestas + Infectadas", "Cuarentena")    

    
    
    figure(9)
    plot(t,matrix(5,:) ,'g')
    title("PERSONAS RECUPERADAS-MUERTAS COVID-19");
    xlabel('$\textbf{Time(days)}$','interpreter','latex')
    ylabel('$\textbf{Number of people}$','interpreter','latex')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin'; 
    hold on
    plot(t,matrix(6,:),'r');
    legend("Recuperadas", "Muertas")    
end
