function [r,jac, ari]=rand_index(U,V,type)
%U,V are membership matrices
NRU=size(U,1);
NRV=size(V,1);
N=size(U,2);
%normalize U,V
% SU=sum(U);
% SV=sum(V);
% U=U./(ones(NRU,1)*SU);
% V=V./(ones(NRV,1)*SV);
 T=U*V';
%TCOLS=sum(sum(T,2));
%TROWS=sum(sum(T,1));
%normalize T
%T=T/TCOLS*TROWS;
    
T1=T.*(T-1);
T2=T.*T;
TCOLS=sum(T,2);
TROWS=sum(T,1);
TSUM=sum(sum(T));
switch type
    case 1
        %rand index
        r=(TSUM*(TSUM-1)+2*sum(sum(T2))-(sum(TCOLS.*TCOLS)+sum(TROWS.*TROWS)))/(TSUM*(TSUM-1));
        jac=0;
        ari=0;
    case 2
        %rand corrected for chance
        A=0.5*(sum(TCOLS.*(TCOLS-1))*sum(TROWS.*(TROWS-1)))/TSUM/(TSUM-1);
        r=(sum(sum(T1))-A)/(0.5*(sum(TCOLS.*(TCOLS-1))+sum(TROWS.*(TROWS-1)))- A);
        jac=0;
        ari=0;
    %%%%%%%%%%%%%%%%%%%%%%%
    case 3
    %Brower Crisp Rand
        %Make the contingency matrix
        U1=U;
        U2=V;
        NMat=U1*U2'; 
        NC1=NRU;
        NC2=NRV;
        n12 = 0;
        for i=1:NC1
            for j=1:NC2
                if(NMat(i,j)>1)
                    n12 = n12 + 0.5*NMat(i,j)*(NMat(i,j)-1);%nchoosek(NMat(i,j),2);
                end
            end
        end
        n1 = 0;
        for i=1:NC1
            T=sum(NMat(i,:));
            if( T > 1 )
                n1 = n1 + 0.5*T*(T-1);%nchoosek(sum(NMat(i,:)),2);
            end
        end
        n2 = 0;
        for i=1:NC2
            T=sum(NMat(:,i));
            if( T > 1 )
                n2 = n2 + 0.5*T*(T-1);%nchoosek(sum(NMat(:,i)),2);
            end
        end
        %Calc a,b,c,d
        a = n12;
        b = n2 - n12;
        c = n1 - n12;
        TT=sum(sum(NMat));
        d= 0.5*TT*(TT-1)- (a+b+c);
        %d = nchoosek(sum(sum(NMat)),2) - (a+b+c);
        r=(a + d) / (a+b+c+d);
        jac = a / (a+b+c);
        ari=(2*(a*d-b*c)) / ( c*c + b*b + 2*a*d + (a+d)*(c+b));
        ari=(ari+1)/2;
        
    %%%%%%%%%%%%%%%%%
    case 4
        %Brower Extension for crisp
        %U1'*U1 gives us the bonding matrix
        U1=U;
        U2=V;
        B1 = U1'*U1;
        B2 = U2'*U2;    
        Bfinal = min(B1,B2');
        fx = sum(sum(Bfinal)) / 2;
        a = fx - length(Bfinal)/2;   
        B1c = 1-B1;   
        Bfinalc = min(B1c,B2');
        b = sum(sum(Bfinalc)) / 2;      
        B2c = 1-B2;    
        Bfinalc = min(B1,B2c);
        c = sum(sum(Bfinalc)) / 2;          
        Bfinalc = min(B1c,B2c);
        d = sum(sum(Bfinalc)) / 2; 
        r=(a + d) / (a+b+c+d);
        jac = a / (a+b+c);
        ari=(2*(a*d-b*c)) / ( c*c + b*b + 2*a*d + (a+d)*(c+b));
        ari=(ari+1)/2;
    %%%%%%%%%%%%%%%%%
    case 5
        U1=U;
        U2=V;
        U1sum = sum(U1);
        U2sum = sum(U2);    
        U1size = size(U1);
        U1n = U1;
        for i=1:U1size(2)
            U1n(:,i) = U1(:,i) / U1sum(i);
        end
        U2size = size(U2);
        U2n = U2;
        for i=1:U2size(2)
            U2n(:,i) = U2(:,i) / U2sum(i);
        end    
        B1 = U1n'*U1n;
        B2 = U2n'*U2n;

        Bfinal = min(B1,B2');
        fx = sum(sum(Bfinal)) / 2;
        a = fx - length(Bfinal)/2;

        B1c = 1-B1;

        Bfinalc = min(B1c,B2');
        b = sum(sum(Bfinalc)) / 2;  

        B2c = 1-B2;

        Bfinalc = min(B1,B2c);
        c = sum(sum(Bfinalc)) / 2;       

        Bfinalc = min(B1c,B2c);
        d = sum(sum(Bfinalc)) / 2; 
        r=(a + d) / (a+b+c+d);
        jac = a / (a+b+c);
        ari=(2*(a*d-b*c)) / ( c*c + b*b + 2*a*d + (a+d)*(c+b));   
        ari=(ari+1)/2;
    case 6
         Umax=max(U);
         for i=1:NRU
             UI{i}=find(U(i,:)==Umax);
         end
         Vmax=max(V);
         for i=1:NRV
             VI{i}=find(V(i,:)==Vmax);
         end
         %build the match matrix (put the max nr of clusters on lines)
         M=min(NRV, NRU);
         match=zeros(NRU, NRV);
         for i=1:NRU
             for j=1:NRV
                 match(i,j)=length(intersect(UI{i},VI{j}));
                 
             end
         end
         [mm,r]=hungarian(N-T);
         r=(M*N-r)/N;
%          %Viterbi- alpha pass
%          alpha=zeros(N1,M);
%          alpha(:,1)=match(:,1);
%          [mm,indel(1)]=max(match(:,1));
%          for i=2:M
%              [maxel, indel(i)]=max(alpha(setdiff([1:N1],indel),i-1));
%              alpha(:,i)=alpha(:,i-1)+maxel;
%          end
%          %beta pass
%          r=max(alpha(:,M))
% %          availableN=[1:N1];
% %          for i=M:-1:1
% %              [rv,ri]=max(alpha(availableN,i));
% %              r=r+rv;
% %              availableN=setdiff(availableN,ri);
% %          end
         jac=0;
         ari=0; 
    otherwise
    disp('bad Rand type!\n');
end
