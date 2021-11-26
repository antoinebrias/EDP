function [M,H]=nph_models(model,p,v)
switch model,%new version with parameters and yield built into handle
    case 1, %staged ricker with DD growth (Neubert&Caswell 2000)
        %p=[Rmax Gmax Sjuv Sadult]
        q1=p(5);
        M=@(x,u) [p(4)*(1-u) p(3)*(1-q1*u)*p(2)*exp(-.01*(x(1)*(1-u)+x(2)*(1-q1*u))); p(1)*(1-u)*exp(v(1)*randn(1)-v(1)^2/2) p(3)*(1-q1*u)*(1-p(2)*exp(-.01*(x(1)*(1-u)+x(2)*(1-q1*u))))]*x;
        H=@(x,u) x(1)*u+x(2)*q1*u;
    case 2, %ricker with two locations
        %p=[Rmax_1 Rmax_2 m1 m2]
        q1=p(5);
        M=@(x,u) [(1-p(3)) p(4);p(3) (1-p(4))]*[p(1)*x(1)*(1-u)*exp(-x(1)*(1-u));p(2)*x(2)*(1-q1*u)*exp(-x(2)*(1-q1*u))].*exp(randn(2,1)*v(1)-v(1)^2/2);
        H=@(x,u) x(1)*u+x(2)*q1*u;
     case 3, %maternal effects sensu Ginzburg
        %p=[Rmax a M]
        M=@(x,u) [p(1)*exp(randn(1)*v(1)-v(1)^2/2)*x(1)*(1-u)*x(2)^p(3)/(1+x(2)^p(3)); 
            p(2)/(1+p(1)*x(1)*(1-u)*x(2)^p(3)/(1+x(2)^p(3)))];
        H=@(x,u) x(1)*u;
    case 4, %'Ricker with seasonal growth'
        %p=[Rmax A 2pi*f]
        M=@(x,u) [x(1)*(1-u)*exp(p(1)*(1+p(2)*x(2))-x(1)*(1-u)+randn(1)*v(1)-v(1)^2/2);sin(asin(x(2))+p(3))];
        H=@(x,u) x(1)*u;
%     case 5, %competition or predation
%         %p=[Rmax_1 Rmax_2 b1 b2];
%         q1=p(5);
%         M=@(x,u) [p(1)*x(1)*(1-u)*exp(-x(1)*(1-u)+p(3)*x(2)*(1-q1*u));p(2)*x(2)*(1-q1*u)*exp(-x(2)*(1-q1*u)+p(4)*x(1)*(1-u))].*exp(randn(2,1)*v(1)-v(1)^2/2);
%         H=@(x,u) x(1)*u+x(2)*q1*u;

    case 5, %Modified Nicholson-Bailey 
        %p=[Rmax_1 Rmax_2 b1 b2];
        q1=p(5);
        M=@(x,u) [p(1)*x(1)*(1-u)*exp(-x(1)*(1-u)-p(3)*x(2)*(1-q1*u));
                  p(2)*x(1)*(1-u)*(1-exp(-p(3)*x(2)*(1-q1*u)))].*exp(randn(2,1)*v(1)-v(1)^2/2);
        H=@(x,u) x(1)*u+x(2)*q1*u;
% %     case 5, %predation
%         %p=[Rmax_1 Rmax_2 b1 b2];
%         q1=p(5);
%          M=@(x,u) [p(1)*x(1)*(1-u)*exp(-x(1)*(1-u)+p(3)*x(2)*(1-q1*u));p(2)*x(2)*(1-q1*u)*exp(+p(4)*x(1)*(1-u))].*exp(randn(2,1)*v(1)-v(1)^2/2);
%         H=@(x,u) x(1)*u+x(2)*q1*u;
    case 6, %contemporary evolution  from Doebeli & deJong TPB 1999 with selectivity
        if length(p)==5,q1=p(4);q2=p(5);else q1=0.8;q2=0.6;end
        nst=8^p(2);
        B=@(n) [n(1)+n(2) n(2)/4 0;
                n(3) n(1)+n(2)/2+n(3) n(1);
                0 n(2)/4 n(2)+n(3)]/(sum(n)+(sum(n)==0));
        S=@(x,p,v) [1/(1+(p(1)-1)/nst*x^p(2));1/(1+(p(1)-1)/nst*(p(3)*x)^p(2));1/(1+(p(1)-1)/nst*x^p(2))].*exp(randn(3,1)*v(1)-v(1)^2/2);
        F=@(x,u) [(1-u); (1-q1*u); (1-q2*u)].*x;
        %fishing first
        M=@(x,u) p(1)*S(sum(F(x,u)),p,v).*(B(F(x,u))*F(x,u));
        H=@(x,u) sum(x-F(x,u));
    case 7, %classic scalar Ricker 
        %p=[Rmax]
        M=@(x,u) p(1)*x*(1-u)*exp(-x*(1-u)+v(1)*randn(1)-v(1)^2/2);
        H=@(x,u) x(1)*u;
    case 8, %delated Hassel
        %p=[Rmax_1 Rmax_2 b1 b2];
        q1=p(5);
        M=@(x,u) [p(1)*x(1)*(1-u)/(1+x(2))^p(2).*exp(randn(1,1)*v(1)-v(1)^2/2);
                  x(1)*(1-u)];
        H=@(x,u) x(1)*u;
    case 9, %delayed logistic
        %p=[Rmax Sjuv Sadult sel_juv]
        q1=p(4);
        M=@(x,u) [min(1,max(0,p(1)*(1-u)*x(1)*(1-x(2))*(1+v(1)*(1-x(2))*(rand(1)-.5)*2))); 
                  (1-u)*x(1)];
        H=@(x,u) x(1)*u;
    case 10, %shepherd competition - this is a simple extension of the shepherd model for 2 species
%         %p=[Rmax_1 Rmax_2 b1 b2];
        q1=p(5);m=1.5;
        M=@(x,u) [p(1)*m*x(1)*(1-u)/(1+(p(1)*m-1)*(x(1)*(1-u)+p(2)*x(2)*(1-q1*u))^p(4));
                  p(1)*x(2)*(1-q1*u)/(1+(p(1)-1)*(p(3)*x(1)*(1-u)+x(2)*(1-q1*u))^p(4))].*exp(randn(2,1)*v(1)-v(1)^2/2);
        H=@(x,u) x(1)*u+x(2)*q1*u;
end


%%old version with parameters as inputs, yield unspecified
% switch model,
%     case 1, %staged ricker with DD growth (Neubert&Caswell 2000)
%         %p=[Rmax Gmax Sjuv Sadult]      
%         M=@(x,p,v,u) [p(4)*(1-u) p(3)*p(2)*exp(-.01*(x(1)*(1-u)+x(2))); p(1)*(1-u)*exp(v(1)*randn(1)-v(1)^2/2) p(3)*(1-p(2)*exp(-.01*(x(1)*(1-u)+x(2))))]*x;
%         Eq=@(p,u) [];
%     case 2, %ricker with two locations
%         %p=[Rmax_1 Rmax_2 m1 m2]
%         M=@(x,p,v,u) [(1-p(3)) p(4);p(3) (1-p(4))]*[p(1)*x(1)*(1-u)*exp(-x(1)*(1-u)+randn(1)*v(1)-v(1)^2/2);p(2)*x(2)*exp(-x(2)+randn(1)*v(1)-v(1)^2/2)];
% %     case 3, %maternal effects sensu Ginzburg
% %         %p=[Rmax a M]
% %         M=@(x,p,v,u) [p(1)*exp(randn(1)*v(1)-v(1)^2/2)*x(1)*(1-u)*x(2)^p(3)/(1+x(2)^p(3)); 
% %             p(2)*x(2)/(1+p(1)*x(1)*(1-u)*x(2)^p(3)/(1+x(2)^p(3)))];
%      case 3, %maternal effects sensu Ginzburg
%         %p=[Rmax a M]
%         M=@(x,p,v,u) [p(1)*exp(randn(1)*v(1)-v(1)^2/2)*x(1)*(1-u)*x(2)^p(3)/(1+x(2)^p(3)); 
%             p(2)/(1+p(1)*x(1)*(1-u)*x(2)^p(3)/(1+x(2)^p(3)))];
%     case 4, %'Ricker with seasonal growth'
%         %p=[Rmax A 2pi*f]
%         M=@(x,p,v,u) [x(1)*(1-u)*exp(p(1)*(1+p(2)*x(2))-x(1)*(1-u)+randn(1)*v(1)-v(1)^2/2);sin(asin(x(2))+p(3))];
%     case 5, %competition or predation
%         %p=[Rmax_1 Rmax_2 b1 b2];
%         M=@(x,p,v,u) [p(1)*x(1)*(1-u)*exp(-x(1)*(1-u)-p(3)*x(2)+randn(1)*v(1)-v(1)^2/2);p(2)*x(2)*exp(-x(2)-p(4)*x(1)*(1-u)+randn(1)*v(1)-v(1)^2/2)];
%         Eq=@(p,u) [log(p(1)*(1-u))-p(3)*log(p(2));log(p(2))-p(4)*log(p(1)*(1-u))]/(1-p(3)*p(4));
% %     case 4, %'Ricker with seasonal K'
% %         %p=[Rmax A 2pi*f]
% %         M=@(x,p,v,u) [p(1)*x(1)*(1-u)*exp(-x(1)*(1-u)*(1+p(2)*x(2))+randn(1)*v(1)-v(1)^2/2);sin(asin(x(2))+p(3))];
% %     case 6, %contemporary evolution  from Doebeli & deJong TPB 1999
% %         M=@(x,p,v,u) [p(1)*x(1)*(1-u)*((x(2)^2+(1-x(2))^2)/(1+(x(1)*(1-u))^p(2))+2*x(2)*(1-x(2))/(1+(p(3)*x(1)*(1-u))^p(2)))*exp(v(1)*randn(1)-v(1)^2/2);
% %             x(2)*(x(2)/(1+(x(1)*(1-u))^p(2))+(1-x(2))/(1+(p(3)*x(1)*(1-u))^p(2)))/((x(2)^2+(1-x(2))^2)/(1+(x(1)*(1-u))^p(2))+2*x(2)*(1-x(2))/(1+(p(3)*x(1)*(1-u))^p(2)))];
%     case 6, %contemporary evolution  from Doebeli & deJong TPB 1999 with selectivity
%         q1=0;q2=0;
%         B=@(n) [n(1)+n(2) n(2)/4 0;
%                 n(3) n(1)+n(2)/2+n(3) n(1);
%                 0 n(2)/4 n(2)+n(3)]/sum(n);
%         S=@(x,p,v) [1/(1+(p(1)-1)*x^p(2));1/(1+(p(1)-1)*(p(3)*x)^p(2));1/(1+(p(1)-1)*x^p(2))].*exp(randn(3,1)*v(1)-v(1)^2/2);
%         F=@(u) [(1-u); (1-q1*u); (1-q2*u)];
%         %fishing first
%         M=@(x,p,v,u) p(1)*S(F(u)'*x,p,v).*(B(F(u).*x)*(F(u).*x));
% end
%         
