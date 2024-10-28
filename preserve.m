function new_P = preserve(LastP,Derta,option)
%{
 Input argument:
        LastP : 最后一次的的控制顶点 2-by-N
        Derta : 最后一次更新控制顶点的Derta  2-by-N
%}
new_P = LastP;
if option
    PerLastP = LastP-Derta; % 倒数第二次的P
    i=2;
    tempA1 = [Derta(:,i),new_P(:,i-1)-PerLastP(:,i+1)];%中间
    tempB1 = new_P(:,i-1)-PerLastP(:,i);
    tempA3 = [Derta(:,i),PerLastP(:,i+2)-PerLastP(:,i+1)];
    tempB3 = PerLastP(:,i+2)-PerLastP(:,i);
    T1 = tempA1\tempB1;
    T3 = tempA3\tempB3;
    T = [T1(1),T3(1)];
    T = min(T(T>=0 & T<=1));
    if T<1e-1
        tT = 0;
    end
    if isempty(T)
        new_P(:,i) = LastP(:,i);
    else
        M = PerLastP(:,i)*(1-T)+T*LastP(:,i);
        new_P(:,i) = PerLastP(:,i)+0.8*(M-PerLastP(:,i));
    end
    for i = 3:size(LastP,2)-2
        tempA1 = [Derta(:,i),new_P(:,i-1)-PerLastP(:,i+1)];%中间
        tempB1 = new_P(:,i-1)-PerLastP(:,i);
        tempA2 = [Derta(:,i),new_P(:,i-1)-new_P(:,i-2)];%左
        tempB2 = new_P(:,i-1)-PerLastP(:,i);
        tempA3 = [Derta(:,i),PerLastP(:,i+2)-PerLastP(:,i+1)];%右
        tempB3 = PerLastP(:,i+2)-PerLastP(:,i);
        
        T1 = tempA1\tempB1;
        T2 = tempA2\tempB2;
        T3 = tempA3\tempB3;
        T = [T1(1),T2(1),T3(1)];
        T = min(T(T>=0 & T<=1));
        if T<0.5
            tT = 0;
        end
        if isempty(T)
            new_P(:,i) = LastP(:,i);
        else
            M = PerLastP(:,i)*(1-T)+T*LastP(:,i);
            new_P(:,i) = PerLastP(:,i)+0.8*(M-PerLastP(:,i));
        end
    end
    i=size(LastP,2)-1;
    tempA1 = [Derta(:,i),new_P(:,i-1)-PerLastP(:,i+1)];%中间
    tempB1 = new_P(:,i-1)-PerLastP(:,i);
    tempA2 = [Derta(:,i),new_P(:,i-1)-new_P(:,i-2)];%左
    tempB2 = new_P(:,i-1)-PerLastP(:,i);
    T1 = tempA1\tempB1;
    T2 = tempA2\tempB2;
    T = [T1(1),T2(1)];
    T = min(T(T>=0 & T<=1));
    if T<1e-1
        tT = 0;
    end
    if isempty(T)
        new_P(:,i) = LastP(:,i);
    else
        M = PerLastP(:,i)*(1-T)+T*LastP(:,i);
        new_P(:,i) = PerLastP(:,i)+0.8*(M-PerLastP(:,i));
    end
end
end