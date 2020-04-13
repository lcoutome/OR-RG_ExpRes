function [T4]=adjs_T(T,t)

% t=t';
% T=T';
% 
% t3=t(1):t(end)+1;
% l = 1;  %t4=0;
% for j=1:length(t)-1
% dif(j) = abs(t(j)-t3(l));
% if dif(j)<=1e-2;
% % I4(l)=I(j);     V4(l)=V(j);     
% t4(l)=t(j); T4(l)=T(j);
% % Adjusting data for Ts
% if (T4(l) <= 1) && (l > 1);
%     T4(l) = T4(l-1);
% end
% l=l+1;
% %disp('ent')
% end
% %if t4(l-1)>=2645
% num2str([j/length(t) t(j) t3(l) dif(j) t4(l-1)],'%2.15f\t')
% % [j T4(l-1)]
% pause
% end
% % I=I4;    V=V4;    
% t2=t4; T2=T4;
% figure
% plot(t2,T2)
% vv=find(T2>=1);
% T2(1:vv(1)-1)=T2(vv(1));

T2 = T(T>1);
t2 = t(T>1);
size(t2)
size(T2)
size(t)
T4 = interp1(t2,T2,t);
% figure(55)
% plot(t,T,t,T4,'linewidth',2)
