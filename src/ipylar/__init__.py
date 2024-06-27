#############################################################
# REQUIRED PACAKGES
#############################################################
from ipylar.version import *
import os
import numpy as np
from scipy.integrate import trapezoid
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from tqdm import tqdm

#############################################################
# INITIALIZATION
#############################################################
print("Loading iPyLAR version",version)

#############################################################
# DATA RETRIEVAL
#############################################################
#Root directory
try:
    FILE=__file__
    ROOTDIR=os.path.abspath(os.path.dirname(FILE))
except:
    FILE=""
    ROOTDIR=os.path.abspath('')

def get_data_path(filename):
    return os.path.join(ROOTDIR, 'data', filename)

def list_data():
    return os.listdir(os.path.join(ROOTDIR, 'data'))

#############################################################
# GLOBAL CONSTANTS
#############################################################
class Const(object):
    # Time-zone
    UTCL=5*3600
    
    # Time constants
    hours=3600
    days=86400
    years=365.25*days

    # Physical constants
    
class Util(object):
    def trapezoid_int(f,a,b,dx):
        xs=np.arange(a,b+dx/2,dx)
        ys=f(xs)
        I=trapezoid(ys,xs)
        return I

    def set_xdates(ax):
        locator=mdates.AutoDateLocator(minticks=5,maxticks=50)
        formatter=mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.tick_params(axis='x',rotation=90)
        ax.margins(0)

    def folding_signal(pars,fint=None, # Mandatory
                    yini=None,yend=None, # Optional
                    plot=False,full_output=False,
                    iper=1,nover=10,ylabel=r"Signal" # Decorative
                    ):
        """
        Creates a folded plot of a time series.

        Parameters:
        pars: Folding parameters:
                t0: initial time [days] (float)
                T: Period [days] (float)

        fint: Interpolating function, fint(t) where t is in seconds, returns the signal at time t.
            It should be an instance of scipy.interpolate.interp1d
            
        iper = 0: number of annual cycles to show (zero, if no annual cycles shown)
        
        nover = 10: Oversampling frequency of interpolation
        
        yini, yend: years of start and end of folding
        
        ylabel = "Signal": ylabel of plot.

        full

        Returns:

        phases: phases (array, Nsam)
        
        fsam: annual cycle average and limits (array, Nsam x 3):
            column 0: 10% quantile
            column 1: median
            column 2: 90% quantile
            
        error: error in folding (signal - fsam)/average
        """

        #Parameters
        t0=pars[0] #days
        T=pars[1] #days

        #Initial time in days
        times=fint.x
        if yini is not None:
            tini=datetime(yini,1,1,0,0,0).timestamp()
        else:
            tini=times[0]
        if yend is not None:
            tend=datetime(yend,1,1,0,0,0).timestamp()
        else:
            tend=times[-1]
        
        #Indexes corresponding to year intervals
        ies=np.arange(len(times))
        mini=ies[times>=tini][0]
        mend=ies[times<=tend][-1]

        #Initial time
        ti=t0+times[mini]/Const.days

        #Interpolating values
        tds=np.linspace(times[mini]/Const.days,times[mend]/Const.days,nover*(mend-mini))
        fints=fint(tds*Const.days)

        #Phases
        phs=np.mod((tds-ti)/T,1)

        #Sampling phases
        Nsam=100
        phsams=np.linspace(0,1,Nsam)
        dph=phsams[1]-phsams[0]

        #Averaging 
        fsam=np.zeros((Nsam,3))
        for i,ph in enumerate(phsams):
            cond=(phs<=(ph+dph))&(phs>(ph-dph))
            fsam[i]=np.quantile(fints[cond],[0.1,0.50,0.90])
        error=abs(fsam[:,2]-fsam[:,0]).sum()/abs(fsam[:,1].mean())

        if plot:
            #Plot folded
            fig=plt.figure()
            ax=fig.gca()
            if iper:
                ax.scatter(phs[::iper],fints[::iper],s=0.1)
            ax.plot(phsams,fsam[:,1],'r-')
            ax.fill_between(phsams,y1=fsam[:,0],y2=fsam[:,2],color='r',alpha=0.3)

            #Decoration
            ax.set_title("Folded signal, $t_0$ = "+f"{t0:.2f} days, T = {T:.2f} days [error {error:.2f}]",fontsize=10)
            ax.set_xlabel(r"Cycle phase, $\phi$")
            ax.set_ylabel(ylabel)
            ax.set_xlim((0,1))
            fig.tight_layout()
            
        if full_output:
            return phsams,fsam,error
        else:
            return error

    def folding_optimal(fint,yini=None,yend=None,nover=30,pdays=30,plot=False):

        #Minimization method
        method="COBYLA"

        #Travel through year finding the best period that minimize variances in folding
        fmin=1e100
        for i,d0 in enumerate(tqdm(np.linspace(0,365,pdays),position=0)):
            sol=minimize(
                Util.folding_signal,
                [d0,365.25],
                method=method,
                args=(fint,yini,yend,False,False,0,nover,""),
                tol=1e-9)
            if sol.fun<fmin:
                fit=sol
                fmin=sol.fun
    
        #Plot best bit
        phsams,fsam,error=Util.folding_signal(fit.x,fint,yini,yend,
                                              plot,True,0,nover,"Folded")
        
        return fit.x,phsams,fsam,error

#############################################################
# LOAD MODULES
#############################################################
from ipylar.basins import *
