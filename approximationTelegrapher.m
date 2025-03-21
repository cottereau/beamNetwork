function [Un,t] = approximationTelegrapher(M,K,T,U0,ind,time,Nt)
% Approximation of the telegrapher's equation with choices of time
% discretization scheme
%
% syntax:  [Un,t] = approximationTelegrapher(M,K,T,U0,ind,time,Nt)
%
% the continuous equation is M*U_tt(t) - K*U(t) = 0, for t in [0,T]
% with vanishing initial conditions and boundary conditions in displacement
%  U(ind,t) = U0(ind,t)
%
% time can be chosen among 'leapfrog', 'central', 'midpoint'
% Nt is the number of time steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretization in space and time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(K,1);
Nt = Nt+1;
t = linspace(0,T,Nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Un = zeros(n,Nt);
Un(ind,1) = U0(t(1));
Vdt = zeros(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time schemes for different cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(time)
    case 'leapfrog'
        for i1 = 2:Nt
            Adt2 = -(M\(K*Un(:,i1-1)));
            Vdt  = Vdt + Adt2;
            Un(:,i1) = Un(:,i1-1) + Vdt;
            Un(ind,i1) = U0(t(i1));
        end
    case 'central'
        A0dt2 = -(M\(K*Un(:,1)));
        for i1 = 2:Nt
            Un(:,i1) = Un(:,i1-1) + Vdt + (A0dt2/2);
            Un(ind,i1) = U0(t(i1));
            A1dt2 = -(M\(K*Un(:,i1)));
            Vdt  = Vdt + (A0dt2+A1dt2)/2;
            A0dt2 = A1dt2;
        end
    case 'midpoint'
        A0dt2 = -(M\(K*Un(:,1)));
        for i1 = 2:Nt
            A1dt2 = -((M+(K/4))\(K*(Un(:,i1-1)+Vdt+(A0dt2/4))));
            Un(:,i1) = Un(:,i1-1) + Vdt + ((A0dt2+A1dt2)/4);
            Un(ind,i1) = U0(t(i1));
            Vdt  = Vdt + (A0dt2+A1dt2)/2;
            A0dt2 = A1dt2;
        end
    otherwise
        error('unknown time scheme')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adding last node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Un = [Un(1,:);Un];
