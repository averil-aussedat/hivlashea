# Paper Malkov-Kudryavtsev
# Numerical computation (RK4) of function a(T)
# which gives an analytical solution for the 1-species case


function rk4solver(T,alpha,rho0,z0,y10,y20)

    ### Numerical parameters
    nbiter=1000000
    dt=T/nbiter

    ### Computation of a(t) thanks to a RK4 method
    # Solving system "dt y(t)=F(t,y(t)) where y=(a, dt a)"
    # Definition of function F for the system
    function F1(a,b)
        return b
    end

    function F2(a,b)
        return alpha*alpha/a^3-1.
    end

    # RK4 scheme
    function RK4_syst(y10,y20,nbiter,dt)
        y1=zeros(nbiter+1)
        y2=zeros(nbiter+1)
        y1[1]=y10
        y2[1]=y20
        for i in 1:nbiter
            k11=F1(y1[i],y2[i])
            k21=F2(y1[i],y2[i])
        
            k12=F1(y1[i]+dt*0.5*k11,y2[i]+dt*0.5*k21)
            k22=F2(y1[i]+dt*0.5*k11,y2[i]+dt*0.5*k21)
        
            k13=F1(y1[i]+dt*0.5*k12,y2[i]+dt*0.5*k22)
            k23=F2(y1[i]+dt*0.5*k12,y2[i]+dt*0.5*k22)
        
            k14=F1(y1[i]+dt*k13,y2[i]+dt*k23)
            k24=F2(y1[i]+dt*k13,y2[i]+dt*k23)

            y1[i+1]=y1[i]+dt/6.0*(k11+2.0*k12+2.0*k13+k14)
            y2[i+1]=y2[i]+dt/6.0*(k21+2.0*k22+2.0*k23+k24)  
        end      
        return y1[nbiter+1],y2[nbiter+1]
    end

    # values of a(T) and dt a(T)
    aT,apT=RK4_syst(y10,y20,nbiter,dt)

    return aT,apT
end
