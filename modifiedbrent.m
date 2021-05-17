function [root,info]=modifiedbrent(func,Int,params) % input @func, Int=[Int.a,Int.b], and params={params.root_tol,params.func_tol,params.maxit}
a=Int.a; 
b=Int.b; 
f_a=func(a); % calculate f(a)
f_b=func(b); % calculate f(b)
f_s=f_b;
if f_a*f_b>=0
    error('Root is not bracketed, may not find solution');
end
if abs(f_a)<abs(f_b) % if |f(a)| < |f(b)| then swap(a,b) and swap(f_a,f_b)
    [a,b]=deal(b,a);
    [f_a,f_b]=deal(f_b,f_a);
end % end if
a_0=a;
b_0=b;
f_b_best=f_b;
c=a;% c := a
f_c=f_a;
s=b_0;
mflag=0; % set mflag
i=0;
if abs(f_b)<=params.func_tol
    root=b;
    info.flag=0;
    return
end
while abs(b-a)>params.root_tol && i<params.maxit% repeat until f(b or s) = 0 or |b − a| is small enough (convergence)
    if f_a~=f_c && f_b~=f_c % if f(a) ≠ f(c) and f(b) ≠ f(c) then
        s=(a*f_b*f_c)/((f_a-f_b)*(f_a-f_c))+(b*f_a*f_c)/((f_b-f_a)*(f_b-f_c))+(c*f_a*f_b)/((f_c-f_a)*(f_c-f_b)); % (inverse quadratic interpolation)
    else % else
        s=b-f_b*(b-a)/(f_b-f_a); % (secant method)
    end % end if
    f_s=func(s);
    if mflag>=5 && abs(b-a)>1/2*abs(b_0-a_0) || abs(f_s)>(1/2*abs(f_b_best))
        s=(a+b)/2;
        mflag=0;
    else
        mflag=mflag+1;
    end
    i=i+1;
    if mflag==0
        f_s=func(s);
    end
    if abs(f_s)<abs(f_b)
        f_b_best=f_s;
    end
    if abs(f_s)<params.func_tol
        root=s;
        info.flag=0;
        return
    end
    %    calculate f(s)
    c=b;%    c := b
    f_c=f_b;
    if (f_a*f_s)<0%    if f(a)f(s) < 0 then
        b=s;%        b := s 
        f_b=f_s;
    else%    else
        a=s;%        a := s
        f_a=f_s;
    end%    end if
    if abs(f_a)<abs(f_b)%    if |f(a)| < |f(b)| then
        [a,b]=deal(b,a);
        [f_a,f_b]=deal(f_b,f_a);     
    end%    end if
    if mflag==0
        a_0=a;
        b_0=b;
    end
end%end repeat
if i>params.maxit
    root=s;
    info.flag=1;
end
root=s;
info.flag=0;
