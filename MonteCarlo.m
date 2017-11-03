eps = 32.4;
sigma = 2.96;
molecules = 2.^(1:9);
T = 300;Kb=1;
rmin = zeros(1,9);Uavg = zeros(1,9);
for N = molecules;
    N
x = (0.25:1.5:19.75)';
A = randi(14,N,3);
A = 1.5*A - 1.25;
Utotal = 0; Utotalcol = zeros(1,9); r1 =Inf;Utotalprime=0;
for i=1:N-1
    for j = i+1:N
        if A(i,:) == A(j,:)
            idz = 1;
            while idz ~= 0               
            np = 1.5*randi(14,1,3)-1.25;
            idx = find( A(:,1) == np(1) );
            idy = find( A(idx,2) == np(2) );
            idz = find( A(idy,3) == np(3) ); 
            end
            A(j,:) = np;
        end 
        r = ((A(i,1)- A(j,1))^2 + (A(i,2)- A(j,2))^2 + (A(i,3)- A(j,3))^2 )^0.5;  
        U = 4*eps*(sigma^12/r^12 - sigma^6/r^6);
        if r<=10
        Utotal = Utotal + U;
        end
        if r<r1
            r1=r;
        end
    end
end
A;
rmin(log(N)/log(2))=r1;
%MC loop
for k = 1:100000
S = randi(N); % select random molecule
SS = 1:N;
SS(S)=[]; %prevent s-s interaction
Uold = 0; Unew = 0;
for i=SS
        r = ((A(i,1)- A(S,1))^2 + (A(i,2)- A(S,2))^2 + (A(i,3)- A(S,3))^2 )^0.5;
        U = 4*eps*(sigma^12/r^12 - sigma^6/r^6);
        if r<=10
        Uold = Uold + U;
        end
end
            newpos = A(S,:) + rand(1,3) - 0.5; %new pos
            while newpos(1)<0 || newpos(2)<0 || newpos(3)<0 || newpos(1)>20 ||newpos(2)>20 ||newpos(3)>20
                 newpos = A(S,:) + rand(1,3) - 0.5;
            end
            
for i=SS
        r = ((A(i,1)- newpos(1))^2 + (A(i,2)- newpos(2))^2 + (A(i,3)- newpos(3))^2 )^0.5;  
        U = 4*eps*(sigma^12/r^12 - sigma^6/r^6);
        if r<=10
        Unew = Unew + U;
        end
end

if rand <= min(1,exp((Uold-Unew)/(T*Kb))) %moving criterion
        A(S,:)=newpos;
    Utotal = Utotal + Unew - Uold ;
    if k>10000
        Utotalprime = Utotalprime + Utotal;
    end
end
if rem(k,10000)== 0
    k
end
end %end of MC loop
Uavg(log(N)/log(2)) = (1/90000)*(Utotalprime);
end