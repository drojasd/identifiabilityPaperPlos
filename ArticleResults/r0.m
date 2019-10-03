function y=r0(r0,eq,level,T,name,xdata,ydata)
    a=symvar(r0);
    len=length(a);    
    for i=1:len
        ce{i}=char(a(i));
    end
    D1 = floor(sqrt(len)); % Number of rows of subplot
    D2 = D1+ceil((len-D1^2)/D1); 
    y=zeros(len,length(xdata));
    figure
    for i=1:len
        assume(a(i)>0)
        temp=subs(r0,a(setdiff(1:len,i)),[eq{ce(setdiff(1:len,i)),name}(:,1)]');
        param=vpa(solve(temp==level,a(i)));
        subplot(D1,D2,i)
        if ~isempty(param)
            T2=T;
            T2{ce(i),'Estlsqc'}(1)=param;
            y(i,:)=gsua_deval(T2{:,'Estlsqc'}(:,1)',T,xdata);
            plot(xdata,y(i,:),xdata,ydata)
        else
            plot(xdata,ydata)
        end
        title(strcat(ce{i},' = ',num2str(round(double(param),2))))
        assume(a(i),'clear')
    end
end
