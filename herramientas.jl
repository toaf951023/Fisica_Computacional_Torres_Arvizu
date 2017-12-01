__precompile__()
module herramientas
export espaaciofaserungeKutta2e
"""Esta funcion sirve para resover númericamente ecuaciones diferenciales por el método de Runge-Kutta, ya sean sistemas o de cualquier orden orden. 
La función recibecomo argumentos (f,t0,t,y0,delta) una función (f),los puntos incial y final del intervalo [t0,tf], las condiciones inciales (y0), y el tamaño de paso (delta). 
La función devuelve un arreglo con el intervalo [t0,tf] , un arreglo con la solución y un arreglo con la derivada en cada punto, es perfecta para hacer espacios fase  """
function espaaciofaserungeKutta2e(f,t0,tf,y0,delta)
   dyarreglo=[]
    yArreglo = []
    tArreglo=[]
    push!(yArreglo,y0)
    push!(tArreglo,t0)
    push!(dyarreglo,f(t0,y0))
    k1(t,y) = f(t,y)
    k2(t,y) = f(t+(delta/2),y+((delta/2)*k1(t,y)))
    k3(t,y) = f(t+(delta/2),y+((delta/2)*k2(t,y)))
    k4(t,y) = f(t+(delta),y+(delta)*k3(t,y))
    h1(t,y)=k1(t,y)+2*k2(t,y)+2*k3(t,y)+k4(t,y)
    while (tArreglo[length(tArreglo)]<tf)
        push!(tArreglo,tArreglo[length(tArreglo)]+delta)
        push!(yArreglo,yArreglo[length(yArreglo)]+(delta/6)*h1(tArreglo[length(tArreglo)-1],yArreglo[length(yArreglo)]))
        push!(dyarreglo,f(tArreglo[length(tArreglo)-1],yArreglo[length(yArreglo)-1]))
    end
   return [tArreglo, yArreglo, dyarreglo]
end

export riemman
""""Esta función da una integral númerica usando la definición integral de Riemman.
Recibe como argumentos (a,b,f,n) el limite inferior (a), el limite superior (b), la función (f),  así como el tamaño de la partición (n)."""


function riemman(a,b, f::Function,n::Int64)

    deltax=(b-a)/n
    S=0.0
    xo=a
    So=0.0
   xn=0.0
    m=1
   for i in 1:n-1
        xm=a+m*deltax
        
        S+=f(xm)
        m=m+1
        xo=xm
        So=S
    
       
    end 
    Sf=S*(b-a)/n
    return Sf
end 
export Simpson
"""Esta función integra numéricamente la función f(x) en el intervalo [a,b] usando la regla  de Simpson. 
Recibe como argumentos (a,b,f,n) el limite inferior (a), el limite superior (b), la función (f),  así como el tamaño de la partición (n). """
function Simpson(a,b, f,n)

    deltax=(b-a)/n
    S=0.0
    xm=0.0
    So=0.0
    x0=a
    x1=0.0
   for i in 1:n-1
        x1=a+i*deltax
        xm=(x1+x0)*(1/2)
        S=(deltax/6)*(f(x0)+4*f(xm)+f(x1))+So
        So=S
        x0=x1
    end
     
    return So
end 
export eulerFunction
""" Esta función aplica el método de Euler independiente de la dimensión a una función.
La función recibe como argumentos (f,time0,timef,y0,delta) la función (f),los puntos incial (time0) y final del intervalo (timef), [t0,tf], las condiciones inciales (y0), y el tamaño de paso (delta). 
La función devuelve un arreglo con el intervalo y un arreglo con la solución, aunque habra que extraer la solución si la dimensión es mayor a 1 .  """
function eulerFunction(f,time0,timef,y0,delta)
    #Se declaran dos arreglos a los cuales iremos agregando los valores calculados de t y y
    yarreglo=[]
    tarreglo=[]
    
    push!(tarreglo,time0)
    push!(yarreglo,y0)
  
    while(tarreglo[length(tarreglo)]<timef)
        
        
        push!(tarreglo,tarreglo[length(tarreglo)]+delta)
        push!(yarreglo,yarreglo[length(yarreglo)]+f(tarreglo[length(tarreglo)-1],yarreglo[length(yarreglo)])*delta)
        
    end
   return [tarreglo,yarreglo]

end
export trapecio
"""Esta función integra numéricamente la función f(x) en el intervalo [a,b] usando la regla  del trapecio. 
Recibe como argumentos (f,a,b,n) que son   el limite inferior (a), el limite superior (b), la función f, así como el tamaño de la partición n."""
function trapecio(a,b, f::Function,n::Int64)

    deltax=(b-a)/n
    S=0.0
    So=f(a)+f(b)
   
   for i in 1:n-1
        xm=a+i*deltax
        
        S+=2*f(xm)
        
        
        So=S
     
       
    end 
    Sf=S*(b-a)/2n
    return Sf
end 
Pkg.add("SymPy")
using SymPy
A,x,a,n,m,xi,x1=symbols("A,x,a,n,m,xi,x1")
export dif
""" Esta función te permite obtener una derivada simbólica respecto a la variable x"""
function dif(f::Function)

 
    df= diff(f(x),x)

return df
end 
export dif2
""" Esta función te permite obtener una derivada simbólica respecto a la variable xpara el método de Euler implítcito"""
function dif2(f::Function)

 
    df= diff(f(x1),x1)

return df
end 
export newton1
""" Esta función aplica el método de Newton si se da una condición inicial."""

function newton1(f,w0) #se crea una funcion que tome los parametros que reuiera la ecuacion par solucionarse  
    w1=w0
    d=herramientas.dif(f)
    df=lambdify( d,[x]);
    w2=0.0
    tolerancia=0.000000000001
    error=0.001
    raices=[]
    while error >  tolerancia
        w2=w1-(f(w1)/df(w1))
        error=abs(w2-w1)
        w1=w2
       
    end 
    return w1
end
export b
""" Esta función retira las raices repetidas en un arreglo"""
function b(x)
n=length(x)
    t=[]
    push!(t,x[1])
    for i in 1:n
       if x[i] ∉ t
          push!(t,x[i])  
    
       end 
    end
   return t 
    
 end
export buscadoraderaices
""" Esta función retira las raices repetidas en un arreglo"""
function buscadoraderaices(x)
n=length(x)
    raicesrepetidas=[]
    y=0.0
    epsilon=0.00000001
    push!(raicesrepetidas,x[1])
    for i in 2:n
         y=x[i]
        v=x[i-1]
          if abs(y-v)>epsilon
         
            push!(raicesrepetidas,y)
        end   
      end 
    return raicesrepetidas
end 
export newton3
"""Método de Newton especial para el método de Euler implícito."""
function newton3(f::Function,w0) #se crea una funcion que tome los parametros que reuiera la ecuacion par solucionarse  
    w1=w0
    d=dif(f)
    df=lambdify( d,[x1]);
    w2=0.0
    tolerancia=0.000000000001
    error=0.001
    raices=[]
    while error >  tolerancia
        w2=w1-(f(w1)/df(w1))
        error=abs(w2-w1)
        w1=w2
       
    end 
    return w1
end
export eulerFunctionimplicita
""" Esta función aplica el método de Euler implícito.
Recibe como argumentos (f,a,b,n) que son la función f, el limite inferior (a), el limite superior (b), así como el tamaño de la partición n.
Devuelve un arreglo con el intervalo y un arreglo con la solución. """
function eulerFunctionimplicita(f,time0,timef,y0,delta)
    #Se declaran dos arreglos a los cuales iremos agregando los valores calculados de t y y
    yarreglo=[]
    tarreglo=[]
    
    push!(tarreglo,time0)
    push!(yarreglo,y0)
  
    while(tarreglo[length(tarreglo)]<timef)
        
        push!(tarreglo,tarreglo[length(tarreglo)]+delta)
        g(x1)=x1-yarreglo[length(yarreglo)]-delta*f(tarreglo[length(tarreglo)],x1)
        y1=newton3(lambdify(g(x1),[x1]),tarreglo[length(tarreglo)])
        push!(yarreglo,y1)
        
    end
    return [tarreglo,yarreglo]
end
export newton2
""" Esta función aplica el método de Newton si se da un arreglo inicial."""
function newton2(f,  w0) #se crea una funcion que tome los parametros que reuiera la ecuacion par solucionarse
   n=length(w0)
   raices=[]
    y=zeros(1,n)
    w1=0.0
      d=herramientas.dif(f)
    t=zeros(1,n)
    df=lambdify( d,[x]);   
    for i in 1:n
        w1=w0[i]

        for j in 1:300
        
             w1=w1-(f(w1)/df(w1))
        
        end 
        push!(raices,w1)
    end 
    y= herramientas.buscadoraderaices(raices)
  t=herramientas.b(y)
    return t
end 

export  newtongneral
""" Esta función aplica el método de Newton si se da un arreglo inicial o punto. Requiere como argumentos (f,w0) que son la función (f) y la condición incial (w0) como entradas (f,w0). En casp"""
function newtongneral(f::Function,  w0)
   solucion=0.0
    if ndims(w0)==0
        solucion=herramientas.newton1(f,w0)
    else
        solucion=herramientas.newton2(f,w0)
    end 
    
    return solucion
end
end
