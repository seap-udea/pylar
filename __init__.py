import numpy as np
from scipy.integrate import odeint,solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rc
rc("animation",html="jshtml")
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings("ignore")
from scipy.interpolate import interp1d
from scipy.integrate import quad 

#############################################
#PHYSICAL ROUTINES
#############################################
def toyLAR(fs,t,*vfuns):
  
  #Storages
  L=fs[0]
  A=fs[1]

  #Valves
  nqfun,Pfun,Efun,Qfun=vfuns
  nq=nqfun(t,L,A)
  P=Pfun(t,L,A)
  E=Efun(t,L,A)
  Q=Qfun(t,L,A)

  #Rates
  dLdt=P-E-Q
  dAdt=-P+E+nq

  return np.array([dLdt,dAdt])

def solveLAR(S0,tini,tend,Nt,vfuns,h=0.01):
  
  #Sampling times
  ts=np.linspace(tini,tend,Nt)
  
  #Solution
  solution=odeint(toyLAR,S0,ts,h0=h,args=vfuns)

  #Extract solutions
  Ls=solution[:,0]
  As=solution[:,1]

  #Functions
  nqfun,Pfun,Efun,Qfun=vfuns

  #Value variables
  nqs=nqfun(ts,Ls,As)
  Qs=Qfun(ts,Ls,As)
  Ps=Pfun(ts,Ls,As)
  Es=Efun(ts,Ls,As)
  
  return ts,Ls,As,nqs,Qs,Ps,Es

def massConservation(solution):
  ts,Ls,As,nqs,Qs,Ps,Es=solution

  #Interpolate solution
  nqfun=interp1d(ts,nqs,kind='cubic')
  Qfun=interp1d(ts,Qs,kind='cubic')

  #Integrate
  Qint,err=quad(Qfun,ts[0],ts[-1],epsrel=1e-7)
  nqint,err=quad(nqfun,ts[0],ts[-1],epsrel=1e-7)

  #Comparte
  input=nqint-Qint
  stored=Ls[-1]+As[-1]
  dmass=stored-input

  return nqint,Qint,input,stored,dmass

#############################################
#GRAPHICAL PARAMETERS
#############################################
#Properties
#Size of the box
width=10
#Vertical offset text
voff=0.2
#Margin boxes 
fmar=1
#Separation box
fsep=2
#Height box
hrep=3
#Margin vertical channels
mver=0.2

#Decorative properties
fs=12
#Fluid in repositories
wc='c'
#Fluid in channels
tc='c'
#Color lines of channels
cc='k'
#Width of walls
lw=2

#Size repositories and channels
uL=hrep
uA=hrep
unq=1
uQ=1
uP=1
uE=1

#############################################
#GRAPHICAL ROUTINES
#############################################
def drawRepo(ax,xy,w,h,s,dir=+1,label=""):
  #Direction
  va="bottom" if dir>0 else "top"
  #Border
  ax.add_patch(Rectangle(xy,w,dir*h,color='none',ec='k',lw=lw))
  #Content
  r=Rectangle(xy,w,dir*s,color=wc,lw=0,zorder=-100)
  ax.add_patch(r)
  #Label
  ax.text((2*xy[0]+w)/2,xy[1]+dir*h+dir*voff,label,
          fontsize=1.5*fs,fontweight='bold',ha='center',va=va)
  return r

def drawChannel(ax,xy,w,s,dir=+1,label=""):
  #Direction
  s*=dir
  va="bottom" if dir<0 else "top"
  #Content
  r=Rectangle(xy,w,s,color=wc,ec=wc,lw=lw,zorder=+100)
  ax.add_patch(r)
  #Borders
  bt,=ax.plot([xy[0],xy[0]+w],[xy[1],xy[1]],'-',
             lw=lw,color=cc,zorder=+100)
  bb,=ax.plot([xy[0],xy[0]+w],[xy[1]+s,xy[1]+s],'-',
             lw=lw,color=cc,zorder=+100)
  #Text
  ax.text(xy[0]+w,xy[1]-dir*voff,label,
          ha="right",va=va,size=fs)
  return r,bt,bb,dir

def drawVChannel(ax,xy,h,s,dir=+1,label=""):
  #Direction
  s*=dir
  h*=dir
  ha="right" if dir>0 else "left"
  va="top" if dir>0 else "bottom"
  #Content
  r=Rectangle(xy,s,h,color=wc,ec=wc,lw=lw,zorder=+100)
  ax.add_patch(r)
  #Borders
  bl,=ax.plot([xy[0],xy[0]],[xy[1],xy[1]+h],'-',
             lw=lw,color=cc,zorder=+100)
  br,=ax.plot([xy[0]+s,xy[0]+s],[xy[1],xy[1]+h],'-',
             lw=lw,color=cc,zorder=+100)
  #Text
  ax.text(xy[0]-dir*voff,xy[1]+h-dir*voff,label,
          ha=ha,va=va,size=fs,zorder=+100)
  return r,bl,br,dir

def plotLevels(L,A,Q,nq,P,E):
  #Figures
  fig=plt.figure(figsize=(6,6))
  ax=fig.gca()

  #Transform to graphical values
  sL=L*uL
  sA=A*uA
  sQ=Q*uQ
  snq=nq*unq
  sP=P*uP
  sE=E*uE

  ##################################################################
  #Land
  ##################################################################
  xyL,wL,hL=[fmar,fmar+hrep],width-4*fmar,hrep
  rL=drawRepo(ax,xyL,wL,hL,sL,dir=-1,label="Land")

  #############
  #Q
  #############
  xyQ,wQ=[xyL[0]+wL,xyL[1]],2*fmar
  rQ,bQ,tQ,dQ=drawChannel(ax,xyQ,wQ,sQ,dir=-1,label=r"Q$\rightarrow$")
  
  #############
  #P
  #############
  xyP,hP=[xyL[0]+mver*wL,xyL[1]],fsep
  rP,lP,iP,dP=drawVChannel(ax,xyP,hP,sP,dir=+1,label=r"P$\downarrow$")

  ##################################################################
  #Atmosphere
  ##################################################################
  xyA,wA,hA=[fmar,fmar+hrep+fsep],width-4*fmar,hrep
  rA=drawRepo(ax,xyA,wA,hA,sA,dir=+1,label="Atmosphere")
  
  #############
  #nabla-q
  #############
  xynq,wnq=[xyA[0]+wA,xyA[1]],2*fmar
  rnq,bnq,tnq,dnq=drawChannel(ax,xynq,wnq,snq,dir=+1,label=r"$-\nabla\cdot q\leftarrow$")
  
  #############
  #E
  #############
  xyE,hE=[xyA[0]+(1-mver)*wA,xyA[1]],fsep
  rE,lE,iE,dE=drawVChannel(ax,xyE,hE,sE,dir=-1,label=r"$\uparrow$E")

  #Decoración
  ax.set_xlim((0,width))
  ax.set_ylim((0,width))
  ax.set_axis_off();
  #ax.grid()
  fig.tight_layout();
  plt.close(fig)

  #Store objects
  fig.lar=((rL,xyL,wL),(rA,xyA,wA),\
         (rQ,xyQ,wQ,bQ,tQ,dQ),\
         (rnq,xynq,wnq,bnq,tnq,dnq),\
         (rP,xyP,hP,iP,lP,dP),\
         (rE,xyE,hE,iE,lE,dE),)

  return fig
         
def updateLevels(fig,L,A,Q,nq,P,E):
  #Get objects
  oL,oA,oQ,onq,oP,oE=fig.lar

  #Transform to graphical values
  sL=L*uL
  sA=A*uA
  sQ=Q*uQ
  snq=nq*unq
  sP=P*uP
  sE=E*uE

  #Update Land
  rL,xyL,wL=oL
  rL.set_height(-sL)

  #Update Atmosphere
  rA,xyA,wA=oA
  rA.set_height(+sA)

  #Update Q
  rQ,xyQ,wQ,bQ,tQ,dQ=oQ
  rQ.set_height(-sQ)
  tQ.set_data([xyQ[0],xyQ[0]+wQ],[xyQ[1]-sQ,xyQ[1]-sQ])

  #Update nq
  rnq,xynq,wnq,bnq,tnq,dnq=onq
  rnq.set_height(+snq)
  tnq.set_data([xynq[0],xynq[0]+wnq],[xynq[1]+snq,xynq[1]+snq])

  #Update P
  rP,xyP,hP,lP,iP,dP=oP
  rP.set_width(sP)
  lP.set_data([xyP[0]+sP,xyP[0]+sP],[xyP[1],xyP[1]+dP*hP])

  #Update E
  rE,xyE,hE,iE,lE,dE=oE
  rE.set_width(-sE)
  iE.set_data([xyE[0]-sE,xyE[0]-sE],[xyE[1],xyE[1]+dE*hE])

def animateLAR(ts,nqs,Ls,As,Qs,Ps,Es,it=-1):
  #Limits
  SLmax=max(Ls)
  SAmax=max(As)
  Qmax=max(Qs)
  nqmax=max(nqs)
  Pmax=max(Ps)
  Emax=max(Es)

  #Levels
  SL=Ls[0]/SLmax
  SA=As[0]/SAmax
  Q=Qs[0]/Qmax
  nq=nqs[0]/nqmax
  P=Ps[0]/Pmax
  E=Es[0]/Emax

  plt.show("off")
  fig=plotLevels(SL,SA,Q,nq,P,E)
  ax=fig.gca()
  tT=ax.text(1,0,"Nt = 0, t = 0",size=14,transform=ax.transAxes,
             ha="right",va="bottom")
  
  def animacion(it):
    SL=Ls[it]/SLmax
    SA=As[it]/SAmax
    Q=Qs[it]/Qmax
    nq=nqs[it]/nqmax
    P=Ps[it]/Pmax
    E=Es[it]/Emax
    #print(f"SL = {SL}, SA = {SA}, Q = {Q}, nq = {nq}, P = {P}, E = {E}")
    updateLevels(fig,SL,SA,Q,nq,P,E);
    tT.set_text(f"Nt = {it}, t = {ts[it]:.3f}")
    return fig.lar

  if it<0:
    plt.close(fig)
    Nt=len(ts)
    anim=animation.FuncAnimation(fig,animacion,frames=Nt,interval=100)
    return anim
  else:
    animacion(it)
    return fig

def plotSolution(ts,nqs,Ls,As,Qs,Ps,Es,ini=0,end=0):

  #Slice vector
  end=len(ts) if end==0 else end

  #Slice
  ts=ts[ini:end]
  nqs=nqs[ini:end]
  Ls=Ls[ini:end]
  As=As[ini:end]
  Qs=Qs[ini:end]
  Ps=Ps[ini:end]
  Es=Es[ini:end]

  #Figure
  fig,axs=plt.subplots(3,2,figsize=(6,6),sharex=True)

  #Plots
  i=0;j=0
  #nablaq
  axs[i,j].plot(ts,nqs)
  axs[i,j].set_ylabel(r"$-\nabla\cdot q$");
  axs[i,j].set_ylim((min(nqs),max(nqs)));
  
  i=0;j=1
  #Q
  axs[i,j].plot(ts,Qs)
  axs[i,j].set_ylabel("$Q$");
  axs[i,j].set_ylim((min(Qs),max(Qs)));
  
  i=1;j=0
  #Land
  axs[i,j].plot(ts,Ls)
  axs[i,j].set_ylabel("SL");
  axs[i,j].set_ylim((min(Ls),max(Ls)));
  
  i=1;j=1
  #Air
  axs[i,j].plot(ts,As)
  axs[i,j].set_ylabel("SA");
  axs[i,j].set_ylim((min(As),max(As)));
  
  i=2;j=0
  #Air
  axs[i,j].plot(ts,Es)
  axs[i,j].set_ylabel("E");
  axs[i,j].set_ylim((min(Es),max(Es)));
  
  i=2;j=1
  #Air
  axs[i,j].plot(ts,Ps)
  axs[i,j].set_ylabel("P");
  axs[i,j].set_ylim((min(Ps),max(Ps)));
  
  for row in axs:
    for ax in row:
      ax.set_xlim((min(ts),max(ts)))

  fig.tight_layout()