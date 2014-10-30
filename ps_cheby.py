import numpy as np
from IPython import display

import matplotlib.pyplot as plt
import matplotlib.backends
plt.switch_backend("TkAgg")

# 1D wave equation with Chebyshev
#$ requires get_cheby_matrix

# parameters
nx=125;     # space dimension
nt=5000;    # time steps
c=20.;      # ac. velocity
dt=.00001;  # time increment
it0=50;     # source time
ncheck=50;  # checkpoint
idisp=25;   # display frequency
ag=.05;     # gaussian in space for source

maxdif=-1.0
x=np.zeros((nx+1), dtype=np.float)
# space coordinates
#for ix=0:nx, x(ix+1)=cos(pi*ix/nx);  end
for ix in range(0, nx+1):
    x[ix]=np.cos(np.pi*ix/nx)
    maxdif=max(x[ix],x[ix-1])

# Stability
#disp(sprintf(' Stability criterium %g (pref. < 0.2) ',c*dt/abs(max(diff(x)))));
print "Stability criterium %g (pref. < 0.2) : ", c*dt/np.fabs(maxdif)

# Frequenz
# source and seismogram
s=np.zeros((nt+1), dtype=np.float)
seis=np.zeros((nt+1), dtype=np.float)
a=input(' Give source frequency (e.g. 1000) (Hz) ')
a=1.0/a
# ricker wavelet period a
t0=4*a/np.pi;

for it in range(0, nt):
    t=it*dt
    v1=np.sqrt(np.pi)/2.0
    v2=(1.0/a*np.pi*(t-t0))
    v3=(1.0/a*np.pi*(t-t0))
    v4=np.exp(-v3*v3)
    s[it-1]=v1*v2*v4
#    s[it]=np.sqrt(np.pi)/2.0*(1.0/a*np.pi*(t-t0))*np.exp(-(1.0/a*np.pi*(t-t0))^2)
s=s/np.max(np.max(s))
#s[nt+1]=0

# initialization of space dependent fields
p=np.zeros((nx+1), dtype=np.float)      # pressure
lp=np.zeros((nx+1), dtype=np.float)     # pressure 1s derivative
lp1=np.zeros((nx+1), dtype=np.float)    # pressure 2nd derivative
pold=np.zeros((nx+1), dtype=np.float)   # p at previous time step
pnew=np.zeros((nx+1), dtype=np.float)   # p at next time step

gauss=np.zeros((nx+1), dtype=np.float)

x0=-.25;
if ag != 0:
    gauss=np.exp(-1/(ag*ag)*(x-x0)*(x-x0))
    gauss=gauss.transpose()

print " Begin time extrapolation ... "

# attempt to calculate derivative with Chebyshef
# initialize coordinates
cx=np.zeros((nx+1), dtype=np.float)
dmx=np.zeros(((nx+1)*(nx+1)))
dmx=np.arange((nx+1)*(nx+1),dtype=np.float).reshape(nx+1,nx+1)

for ix in xrange(0, nx):
    x[ix]=np.cos(np.pi*ix/nx)
#print x

# initialize derivative operators  (nx)
cx[0]=2.0;
cx[nx]=2.0;

for i in xrange(1, nx):
    cx[i]=1.0

# diagonal
for i in range(0, nx+1):
    for j in range(0, nx+1):
        if i==j:
            if i!=0:
                if i!=nx:
                    dmx[i,i]=-x[i]/(2.0*(1.0-x[i]*x[i]));
        else:
            dmx[i,j]=(cx[i]*np.power(-1,i+j))/(cx[j]*(x[i]-x[j]));
    
    #  corners
    dmx[0,0]=(2*nx*nx+1)/6;
    dmx[nx,nx]=-dmx[0,0];

fig=plt.figure()

plt.subplot(2,2,1)
plt.plot(np.arange(0,nt+1,1)/2,s)
plt.title('Source time function')

plt.subplot(2,2,2)
plt.plot(x,gauss)
plt.title(' Gauss in space')

plt.subplot(2,1,2)
plt.ylim(-a/10,a/10)
line,=plt.plot(x,p)
plt.xlabel('Distance (km)')

fig.show()

Ddot=np.dot(dmx,dmx)

for it in xrange(0, nt):
    # loop over time order
    # 2nd derivative (D twice operated)
    
    #lp=D*D*p.transpose()
    lp=np.dot(Ddot,np.transpose(p))
    
    if np.isnan(lp).any():
        print "NAN in lp at it = ", it
        exit()
    
	# add source
    lp=lp+gauss*s[it]
    if np.isnan(lp).any():
        print "NAN in lp at it = ", it
        exit()
    
    pnew=(2*p-pold)+(c*c*dt*dt*np.transpose(lp))
    if np.isnan(pnew).any():
        print "NAN in pnew at it = ", it
        exit()
    
    pold=p
    p=pnew
    # set boundaries 0
    p[0]=0
    p[nx]=0

    if np.mod(it,10)==0:
        line.set_ydata(p)
        plt.title('Time step : %g ' %it)
        fig.canvas.draw()python

