#############################################################
# REQUIRED PACAKGES
#############################################################
from ipylar import *
import pandas as pd
import numpy as np
from copy import deepcopy
from tqdm import tqdm
from datetime import datetime
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from dateutil.relativedelta import relativedelta
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter

class PlotProperties(object):
    pass

#############################################################
# CLASS DEVELOPMENT
#############################################################
class Basin(object):
    def __init__(self, 
                 key='basin',
                 name='Basin'):
        
        # Basic properties
        self.key = key
        self.name = name

        # Mandatory columns for basin data
        self.data = pd.DataFrame(
            columns=[
                'datestd', 'datetime', 'time', 
                'R', 'dR', # Discharge and absolute error 
                'nq', 'dnq', # Moisture convergence and absolute error
                'Rref', 'nqref', # Reference discharge and moisture convergence
                ]
        )

        # Set plotProperties
        pp = PlotProperties()
        
        #Safe colors
        pp.Qcolor="#d95f02"
        pp.Rcolor="#7570b3"

        #Common strings
        pp.nam_nablaq="Moisture convergence"
        pp.str_nablaq=r"$Q$"
        pp.nam_R="River discharge"
        pp.str_R=r"$R$"
        pp.sty_nq=dict(color=pp.Rcolor,alpha=1)
        pp.sty_R=dict(color=pp.Qcolor,alpha=1)
        pp.lab_nq=rf"{pp.nam_nablaq} ({pp.str_nablaq})"
        pp.lab_R=rf"{pp.nam_R} ({pp.str_R})"

        #Units
        pp.nam_flu="Flux"
        pp.UQ=1e9/(365.25*86400) #dam^3 -> m^3
        pp.uflux=r"$10^3$ km$^3$/yr"
        #Time
        pp.utime="d"
        pp.UT=1
        #Total stored
        pp.UI=1e3 #km^3 -> 10^3 km^3

        # Accumulation, release
        pp.acc_color='#66c2a5'
        pp.rel_color='#fc8d62'

        self.pp = pp


    def read_basin_data(self,source='era5',rsource=None):

        self.qsource = source
        self.rsource = source

        #Read Data
        filepath=get_data_path(f"{self.key}_RQ_{source}.csv")
        data=pd.read_csv(filepath)
        
        #Dates
        data["datestd"]=data["date"].apply(lambda d:datetime.strptime(d,"%Y-%m-%d").strftime("%Y-%m-16"))
        data["datetime"]=data["datestd"].apply(lambda d:datetime.strptime(d,"%Y-%m-%d"))
        data["time"]=data["datetime"].apply(lambda d:d.timestamp()+Const.UTCL)
        
        #Numeric columms
        columns=['water_vap','water_liq','water_ice']
        data.dropna(subset=columns,inplace=True)
        try:
            for col in columns:
                data[col]=data[col].str.replace(",",".").astype("float")
        except:
            pass
        for col in columns:
            data[col]=pd.to_numeric(data[col],errors='coerce').interpolate()

        # Atmospheric data
        data_atm=data.copy()
        data_atm.set_index(data["datestd"],inplace=True)
        data_atm["nq_mass"]=data_atm[columns].sum(axis=1)
        data_atm["nq"]=data_atm["nq_mass"]/1e3
        
        #Vanilla error
        self.vnq=0.05 
        data_atm["dnq"]=self.vnq*abs(data_atm["nq_mass"])/1e3

        # River data
        if rsource is None:
            columns=["discharge"]
        else:
            source = rsource[0]
            self.rsource = source
            format = rsource[1]
            column = rsource[2]
            filepath=get_data_path(f"{self.key}_R_{source}.{format}")
            if format=="csv":
                data=pd.read_csv(filepath)
            elif format=="xlsx":
                data=pd.read_excel(filepath)
            else:
                raise ValueError("Format not supported")
            data["datestd"]=data.apply(lambda x:f"{x.V1:.0f}-{x.V2:02.0f}-16",axis=1)
            columns=[column]

        #Date
        data["datetime"]=data["datestd"].apply(lambda d:datetime.strptime(d,"%Y-%m-%d"))
        data["time"]=data["datetime"].apply(lambda d:d.timestamp()+Const.UTCL)
        
        #Numeric columns
        data.dropna(subset=columns,inplace=True)
        try:
            for col in columns:
                data[col]=data[col].str.replace(",",".").astype("float")
        except:
            pass
        for col in columns:
            data[col]=pd.to_numeric(data[col],errors='coerce').interpolate()
        
        #Consolidate data
        data_riv=data.copy()
        data_riv.set_index(data["datestd"],inplace=True)
        data_riv["R"]=data_riv[columns].sum(axis=1)
        
        #Vanilla error
        self.vR=0.15 
        data_riv["dR"]=self.vR*abs(data_riv["R"])

        # Merge data
        data_merge=pd.DataFrame()
        merge=data_atm.merge(data_riv,left_index=True,right_index=True)
        data_merge["datestd"]=merge["datestd_x"]
        data_merge.set_index(data_merge["datestd"],inplace=True)
        data_merge["datetime"]=merge["datetime_x"]
        data_merge["time"]=merge["time_x"]
        data_merge["R"]=merge["R"]
        data_merge["dR"]=merge["dR"]
        data_merge["nq"]=merge["nq"]
        data_merge["dnq"]=merge["dnq"]

        # Reference
        data_merge["Rref"]=merge["R"]
        data_merge["nqref"]=merge["nq"]
        
        # Save data
        self.data_merge=data_merge
        self.data=data_merge
        self.data.set_index("datestd",drop=False,inplace=True)

        # Basin properties
        self._calc_basin_properties()

        # Interpolation
        self._calc_basin_interpolation()

        print(self)

    def _calc_basin_properties(self):

        # Time properties
        self.tini=float(self.data_merge["time"].iloc[0])
        self.tend=float(self.data_merge["time"].iloc[-1])
        self.dini=self.data_merge["datetime"].iloc[0]
        self.tspan=self.data_merge["time"].iloc[-1]-self.data_merge["time"].iloc[0]
        self.yspan=self.tspan/Const.years
        self.dspan=self.tspan/Const.days

        # Other properties
        self.ndata=len(self.data_merge)
        self.nmerge=len(self.data_merge)
        
        #Covariance matrix and correlation coefficient
        """
        Cov = 
        [ 
        sigma_x^2  rho_xy sigma_x sigma _y
        rho_xy sigma_x sigma _y  sigma_y^2
        ]
        """
        Cov=np.cov(self.data["R"],self.data["nq"])
        """
        rho_xy = Cov[0,1]/(sigma_x*sigma_y)
        sigma_x = Cov[0,0]**0.5
        sigma_y = Cov[1,1]**0.5
        """
        self.rhoRnq=Cov[0,1]/(Cov[0,0]*Cov[1,1])**0.5
        
        #Random seed
        self.seed=hash(self.key)%(2**32-1)
            
    def set_vanilla_errors(self,vQ=None,vR=None):
        self.vQ=vQ
        self.vR=vR
        if vQ:
            self.data["dnq"]=vQ*abs(self.data["nq"])
        if vR:
            self.data["dR"]=vR*abs(self.data["R"])


    def _calc_basin_interpolation(self):
        #Interpolation
        self.Rint=interp1d(self.data["time"],self.data["R"],kind='slinear')
        self.nqint=interp1d(self.data["time"],self.data["nq"],kind='slinear')

    def __str__(self):
        string = ''
        string += '\n' + f"{self.name} (key: {self.key})"
        string += '\n' + f"Number of data points: merge = {self.nmerge}, total = {self.ndata}"
        string += '\n' + f"Initial time: {self.tini}"
        string += '\n' + f"Initial date: {self.dini}"
        string += '\n' + f"Total data points: {self.ndata}"
        string += '\n' + f"Time span: {self.yspan} a = {self.dspan} d"
        string += '\n' + f"Signal correlation (R - nq): {self.rhoRnq}"
        return string
    
    def __repr__(self):
        return self.__str__()
        
    def shuffle_basin(self):
        """
        Errors assuming normal distribution
        """
        shuffled=deepcopy(self)
        
        #Covariance matrix
        for index in self.data.index:
            row=self.data.loc[index]
            rho=self.rhoRnq
            cov=np.array([
                [row.dR**2,rho*row.dR*row.dnq],
                [rho*row.dR*row.dnq,row.dnq**2]
            ])
            ran=np.random.multivariate_normal([row.R,row.nq],cov,1)
            shuffled.data.loc[index,"R"]=ran[0][0]
            shuffled.data.loc[index,"nq"]=ran[0][1]

        # Derivative properties of the basin
        shuffled._calc_basin_properties()
        shuffled._calc_basin_interpolation()

        return shuffled

    def plot_basin_series(
            self,
            ti_win=2000,dt_win=5,
            over_interp=10,
            show_ref=False):

        #Full signal
        tfulls=np.linspace(self.data["time"].iloc[0],self.data["time"].iloc[-1],over_interp*self.ndata)
        dfulls=[datetime.fromtimestamp(t) for t in tfulls]
        Rfulls=self.Rint(tfulls)
        nqfulls=self.nqint(tfulls)

        #Plot full signal
        fig,axs=plt.subplots(2,1,figsize=(9,8))

        #Full panel
        ax=axs[0]

        ax.plot(self.data["datetime"],self.data["nq"]/self.pp.UQ/self.pp.UI,label=self.pp.lab_nq,**self.pp.sty_nq)
        ax.plot(self.data["datetime"],self.data["R"]/self.pp.UQ/self.pp.UI,label=self.pp.lab_R,**self.pp.sty_R)

        if show_ref:
            ax.plot(self.data["datetime"],self.data["nqref"]/self.pp.UQ/self.pp.UI,label=f"{self.pp.str_nablaq} ref.",**self.pp.sty_nq,ls=':')
            ax.plot(self.data["datetime"],self.data["R"]/self.pp.UQ/self.pp.UI,label=f"{self.pp.str_R} ref.",**self.pp.sty_R,ls=":")
            
        ax.set_ylabel(f"Flux [{self.pp.uflux}]");
        legend=ax.legend(loc=(0.0,1),ncol=6,frameon=False);

        #Zoom
        ax=axs[1]
        ax.plot(dfulls,nqfulls/self.pp.UQ/self.pp.UI,label=rf"{self.pp.nam_R} ({self.pp.str_R})",**self.pp.sty_nq)
        ax.errorbar(self.data["datetime"],self.data["nq"]/self.pp.UQ/self.pp.UI,yerr=self.data["dnq"]/self.pp.UQ/self.pp.UI,ls='None',capsize=5,**self.pp.sty_nq)
        ax.plot(dfulls,Rfulls/self.pp.UQ/self.pp.UI,label=rf"{self.pp.nam_R} ({self.pp.str_R})",**self.pp.sty_R)
        ax.errorbar(self.data["datetime"],self.data["R"]/self.pp.UQ/self.pp.UI,yerr=self.data["dR"]/self.pp.UQ/self.pp.UI,ls='None',capsize=5,**self.pp.sty_R)
        
        if show_ref:
            #Interpolation
            self.Rint=interp1d(self.data["time"],self.data["Rref"],kind='slinear')
            self.nqint=interp1d(self.data["time"],self.data["nqref"],kind='slinear')

            #Full signal
            tfulls=np.linspace(self.data["time"].iloc[0],self.data["time"].iloc[-1],over_interp*self.ndata)
            dfulls=[datetime.fromtimestamp(t) for t in tfulls]
            Rfulls=self.Rint(tfulls)
            nqfulls=self.nqint(tfulls)

            ax.plot(dfulls,nqfulls/self.pp.UQ/self.pp.UI,**self.pp.sty_nq,ls=":")
            ax.errorbar(self.data["datetime"],self.data["nqref"]/self.pp.UQ/self.pp.UI,yerr=self.data["dnq"]/self.pp.UQ/self.pp.UI,ls='None',capsize=5,**self.pp.sty_nq)
            ax.plot(dfulls,Rfulls/self.pp.UQ/self.pp.UI,**self.pp.sty_R,ls=":")
            ax.errorbar(self.data["datetime"],self.data["Rref"]/self.pp.UQ/self.pp.UI,yerr=self.data["dR"]/self.pp.UQ/self.pp.UI,ls='None',capsize=5,**self.pp.sty_R)

        ax.set_ylabel(f"Flux [{self.pp.uflux}]");
        ti_win=datetime(ti_win,1,1)
        dt_win=relativedelta(years=dt_win)
        ax.set_xlim(ti_win,ti_win+dt_win)

        #Source
        ax.text(0.98,0.98,f"Source: {self.rsource}/{self.qsource}",
                ha='right',va='top',transform=ax.transAxes,fontsize=8)

        #Decoration
        for ax in axs:
            locator=mdates.AutoDateLocator(minticks=5,maxticks=50)
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)
            ax.tick_params(axis='x', rotation=90)
            ax.margins(0)

        axs[0].set_title(f"{self.name}",pad=30)
        fig.tight_layout()  

        return fig

    def detect_accumulation_release(self,istart=1,verbose=False,advance=False):    

        ###############################
        # CROSSINGS
        ###############################
        # Storage function
        def AlmFun(t):
            A=self.nqint(t)-self.Rint(t)
            return A

        # Intersection points
        ts=np.array(self.data["time"])
        Qs=np.array(self.data["nq"])
        Rs=np.array(self.data["R"])
        ds=Qs-Rs

        sgnp=np.sign(ds[0])
        trs=[-1e100]
        for i,d in enumerate(ds):
            sgn=np.sign(d)
            if (sgn*sgnp<0):
                tr=-ds[i-1]*(ts[i]-ts[i-1])/(ds[i]-ds[i-1])+ts[i-1]
                if (tr-trs[-1])>50*Const.days:
                    trs+=[tr]
            i+=1
            sgnp=sgn
        trs=np.array(trs[1:])

        # Selected points
        trs=trs[istart:]
        NR=len(trs)
        if verbose:
            print(f"Found {NR} intersections")
        self.NR=NR
        
        ###############################
        # COMPUTE ACCUMULATION
        ###############################
        # Evaluate the type of period
        sgnp=0
        tacus=[]
        Iacus=[]
        trels=[]
        Irels=[]
        if verbose or not advance:
            ies=range(1,NR)
        else:
            ies=tqdm(range(1,NR))

        for i in ies:
            evalue=False

            #Compute accumulation
            I=Util.trapezoid_int(AlmFun,trs[i-1],trs[i],Const.days)/1e9 #10^9 = (10^3)^3 (km^3) 

            #Determine the sign of the integral
            sgn=np.sign(I)
            if i==1:
                if verbose:
                    print("First value")
                sgnp=sgn
                evalue=1
            #Change of sign respect to previous
            elif sgn*sgnp<0:
                evalue=1
                if verbose:
                    print(f"Change of sign in {i} with integral {I}")
            #No change of sign
            else:
                evalue=2
                if verbose:
                    print(f"No change of sign in {i} with integral {I}")

            #If accumulate
            if evalue:
                if sgn>0:
                    #If evalue = 2 is because the previous has the same sign
                    if evalue==1:
                        tacus+=[trs[i-1]]
                        Iacus+=[I]
                    else:
                        Iacus[-1]+=I
                    if verbose:
                        print(f"Storing acummulation integral {i} (t in [{trs[i-1]},{trs[i]}]): {Iacus[-1]}")
                else:
                    #If evalue = 2 is because the previous has the same sign
                    if evalue==1:
                        trels+=[trs[i-1]]
                        Irels+=[I]
                    else:
                        Irels[-1]+=I
                    if verbose:
                        print(f"Storing release integral {i} (t in [{trs[i-1]},{trs[i]}]): {Irels[-1]}")
            sgnp=sgn

        #Check the latest point
        I=Util.trapezoid_int(AlmFun,trs[-2],trs[-1],Const.days)/1e9
        if sgn>0:
            trels+=[trs[-1]]
            Irels+=[I]
        else:
            tacus+=[trs[-1]]
            Iacus+=[I]

        trels=np.array(trels)
        tacus=np.array(tacus)
        Irels=np.array(Irels)
        
        #Restrict trels to intervals of tacus
        cond=(trels>tacus[0])&(trels<tacus[-1])
        trels=trels[cond]
        Irels=Irels[cond]
        
        ###############################
        # COMPUTE ACCUMULATION
        ###############################
        S=0
        Ss=[]
        tas=[]
        dtas=[]
        trs=[]
        dtrs=[]
        tns=[]
        for i in range(len(tacus)-1):
            #Accumulation
            t=(tacus[i]+trels[i])/2
            tas+=[t]
            dta=abs(trels[i]-tacus[i])/86400
            dtas+=[dta]
            Ia=Iacus[i]
            if verbose:
                print(f"Accumulation {i} (t in [{tacus[i]},{trels[i]}]) : {Iacus[i]}")
            #Release
            t=(trels[i]+tacus[i+1])/2
            trs+=[t]
            dtr=abs(tacus[i+1]-trels[i])/86400
            dtrs+=[dtr]
            Ir=Irels[i]
            if verbose:
                print(f"Release {i} (t in [{trels[i]},{tacus[i+1]}]) : {Irels[i]}")
            #Net
            t=(tacus[i]+tacus[i+1])/2
            tns+=[t]
            S+=Ia+Ir
            Ss+=[S]
        
        self.nper=i+1
        self.tns=np.array(tns)
        self.dtns=[datetime.fromtimestamp(t) for t in self.tns]
        self.Ss=np.array(Ss)
        
        #Auxiliar
        self.trs=np.array(trs)
        self.tas=np.array(tas)
        self.dtas=np.array(dtas)
        self.dtrs=np.array(dtrs)
        self.Iacus=np.array(Iacus)
        self.Irels=np.array(Irels)
        self.trels=trels
        self.tacus=tacus

        self.detection=True

    def plot_accumulation_release(self, color='r', drgs=[], basinsets=[]):
        
        if not self.detection:
            print("You cannot plot accumulation/release/storage of a basin whose periods has not been detected. Use detect_accumulation_release(basin).")
            return
        
        #Basin data
        Qs=np.array(self.data["nq"])
        Rs=np.array(self.data["R"])
        
        ##########################################################
        # PERIODS
        ##########################################################
        #Plot
        Np=2
        fig,axs=plt.subplots(Np,1,figsize=(8,2.5*Np))

        dt=self.tspan/Np
        for i,ax in enumerate(axs):
            ti=self.tini+i*dt
            cond=(self.data["datetime"]>=datetime.fromtimestamp(ti))&(self.data["datetime"]<datetime.fromtimestamp(ti+dt))
            ax.plot(self.data[cond]["datetime"],Rs[cond]/self.pp.UQ/self.pp.UI,'-',label=self.pp.lab_R,**self.pp.sty_R)
            ax.plot(self.data[cond]["datetime"],Qs[cond]/self.pp.UQ/self.pp.UI,'-',label=self.pp.lab_nq,**self.pp.sty_nq)
            ax.set_ylabel(f"""Exchange
    ({self.pp.uflux})""",fontsize=8)
            Util.set_xdates(ax)
            for ta,tr in zip(self.tacus[:-1],self.trels):
                ax.axvspan(datetime.fromtimestamp(ta),datetime.fromtimestamp(tr),color=self.pp.acc_color,alpha=0.3)
            for tr,ta in zip(self.trels,self.tacus[1:]):
                ax.axvspan(datetime.fromtimestamp(tr),datetime.fromtimestamp(ta),color=self.pp.rel_color,alpha=0.3)
            ax.set_xlim(datetime.fromtimestamp(ti),datetime.fromtimestamp(ti+dt))
            ax.tick_params(axis='x',labelsize=8)
            ax.margins(0)

        #Decoration
        axs[0].set_title(f"{self.name}",pad=30,fontsize=10)
        axs[0].legend(loc=(0.0,1),ncol=5,frameon=False,fontsize=8)
        
        fig.tight_layout()
        fig_acrel = fig
        
        ##########################################################
        # STORAGE
        ##########################################################
        fig,axs=plt.subplots(3,1,figsize=(8,6),sharex=True)

        ax=axs[0]
        ax.patch.set_visible(False)
        for tr in self.trels:
            ax.axvline(datetime.fromtimestamp(tr),color='k',lw=0.5,alpha=0.3)
        for ta in self.tacus:
            ax.axvline(datetime.fromtimestamp(ta),color='k',lw=0.5,alpha=0.3)

        ax.plot(self.data["datetime"],Rs/self.pp.UQ/self.pp.UI,label=self.pp.lab_R,**self.pp.sty_R,zorder=100)
        ax.plot(self.data["datetime"],Qs/self.pp.UQ/self.pp.UI,label=self.pp.lab_nq,**self.pp.sty_nq,zorder=100)
        
        fmin,fmax=ax.get_ylim()
        fmean=np.mean([(Rs/self.pp.UQ/self.pp.UI).mean()])
        df=max(abs(fmax-fmean),abs(fmin-fmean))
        ax.set_ylim(fmean-df,fmean+df)
        
        ax.legend(loc=(0.0,1),ncol=5,frameon=False);
        ax.set_ylabel(f"""Exchange
    ({self.pp.uflux})""",fontsize=8)
        ax.margins(0)

        axt=ax.twinx()
        axt.set_zorder(-1)
        ax=axs[1]
        axa=axs[2]
        for i in range(len(self.tacus)-1):
            #Accumulation
            axt.bar(datetime.fromtimestamp(self.tas[i]),
                    self.Iacus[i]/self.pp.UI,width=self.dtas[i],color=self.pp.acc_color,alpha=0.5)
            #Release
            axt.bar(datetime.fromtimestamp(self.trs[i]),
                    self.Irels[i]/self.pp.UI,width=self.dtrs[i],color=self.pp.rel_color,alpha=0.5)
            #Net
            ax.bar(datetime.fromtimestamp(self.tns[i]),
                (self.Iacus[i]+self.Irels[i])/self.pp.UI,
                width=self.dtas[i]+self.dtrs[i],color=color,alpha=1)

        #Acummulated
        lss=[":","--","-."]
        if len(basinsets)>0:
            nbsets=len(basinsets)
            for i,basinset in enumerate(basinsets):
                axa.fill_between(basinset.dreals,
                                basinset.Smins/self.pp.UI,basinset.Smaxs/self.pp.UI,
                                color=basinset.color,alpha=0.2,
                                ls=lss[i])
                axa.fill_between(basinset.dreals,basinset.Smins/self.pp.UI,basinset.Smins/self.pp.UI,
                        color=basinset.color,ls=lss[i],
                        alpha=0.2*(nbsets-i),label=f"Relative error in $R$, {100*basinset.vR:.0f}%")
                
            axa.legend(fontsize=8)
            
        axa.plot([datetime.fromtimestamp(t) for t in self.tns],self.Ss/self.pp.UI,'-',color=color)

        ax.set_ylabel(rf"""Rate of storage change
    ($10^3$ km$^3$/yr)""",fontsize=8)
        axt.set_ylabel(r"""Storage change
    ($10^3$ km$^3$)""",fontsize=8)
        axa.set_ylabel(rf"""Accumulated storage 
    ($10^3$ km$^3$)""",fontsize=8)
        axt.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axt.tick_params(axis='y',which='major',labelsize=8)
        #axt.set_zorder(-1)
        Imax,Imin=max(max(self.Iacus)/self.pp.UI,max(self.Irels)/self.pp.UI),min(min(self.Iacus)/self.pp.UI,min(self.Irels)/self.pp.UI)
        dI=max(abs(Imax),abs(Imin))
        axt.set_ylim(-dI,+dI)

        #Source
        """
        ax.text(0.98,0.98,f"Source: {self.rsource}/{self.qsource}",
                ha='right',va='top',transform=ax.transAxes,fontsize=8)
        """
        ax=axs[0]
        ax.bar(datetime.fromtimestamp(self.tacus[0]),self.Iacus[0]/self.pp.UI,width=0,
            color=self.pp.acc_color,alpha=1,label="Accumulation")
        ax.bar(datetime.fromtimestamp(self.trels[0]),self.Irels[0]/self.pp.UI,width=0,
            color=self.pp.rel_color,alpha=1,label="Release")

        labelx=-0.06
        labels="abcdef"
        for i,ax in enumerate(list(axs)):
            ax.yaxis.set_label_coords(labelx, 0.5)
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.tick_params(axis='y',which='major',labelsize=8)
            ax.text(1,1.02,labels[i],
                    fontsize=12,ha='right',va='bottom',transform=ax.transAxes,fontweight="bold")
        
        #DROUGHTS
        for y in drgs:
            axs[1].text(datetime(y,1,1,0,0,0),0.05,f"{y} Drought",color='k',zorder=100,rotation=90,fontsize=8)
        
        #Limits
        ax.set_xlim(datetime.fromtimestamp(self.tacus[0]),
                    datetime.fromtimestamp(self.tacus[-1]))
        ax.tick_params(axis='x',labelsize=8)

        #Decoration
        axs[0].legend(loc=(0.0,1),ncol=5,frameon=False,fontsize=8);
        axs[2].text(0.5,0.1,self.name,fontsize=12,transform=ax.transAxes,
                va='bottom',ha='center')
        
        ax=axs[2]
        locator=mdates.AutoDateLocator(minticks=5,maxticks=50)
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.tick_params(axis='x', rotation=90)
        
        for ax in axs[1],axs[2]:
            ax.margins(0.02)  
        
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0,hspace=0.2)
        fig_stor = fig

        return fig_acrel, fig_stor

    def get_annual_cycle(self,pdays=5,plot=True):
        
        # Get data
        ts=np.array(self.data["time"])
        Qs=np.array(self.data["nq"])/self.pp.UQ/self.pp.UI
        Rs=np.array(self.data["R"])/self.pp.UQ/self.pp.UI
        self.Qmean=Qs.mean()
        self.Rmean=Rs.mean()

        kind_interp="slinear"
        Rint=interp1d(ts,Rs,kind=kind_interp)
        nqint=interp1d(ts,Qs,kind=kind_interp)

        # Perform optimal folding
        self.fit_nq,phsams_nq,fsam_nq,error_nq=Util.folding_optimal(nqint,pdays=pdays)
        self.fit_R,phsams_R,fsam_R,error_R=Util.folding_optimal(Rint,pdays=pdays)

        # Get results for folding
        self.phsams_nq,self.fsam_nq,self.error_nq = Util.folding_signal(
            [0,self.fit_nq[1]],fint=nqint,plot=plot,full_output=True,ylabel='Moisture Convergence'
            )

        self.phsams_R,self.fsam_R,self.error_R = Util.folding_signal(
            [0,self.fit_R[1]],fint=Rint,plot=plot,full_output=True,ylabel='Discharge'
            )

        # Report
        print(f"Average discharge (<R>): {self.Rmean:.2f}")
        print(f"Average moisture convergence (<Q>): {self.Qmean:.2f}")
        print(f"R optimal folding: t = {self.fit_R[0]:.2f} days, P = {self.fit_R[1]:.3f} days, error = {error_R:.2f}")
        print(f"Q optimal folding: t = {self.fit_nq[0]:.2f} days, P = {self.fit_nq[1]:.3f} days, error = {error_nq:.2f}")

    def plot_annual_cycle(self):

        if 'fit_R' not in self.__dict__.keys():
            raise ValueError("You must get annual cycle first")

        #############################################
        # ANNUAL CYCLE PLOT
        #############################################
        fig,axs=plt.subplots(1,1)
        ax=axs
        ax.plot(self.phsams_R,self.fsam_R[:,1],color=self.pp.Rcolor,label="River discharge (R)")
        ax.fill_between(self.phsams_R,self.fsam_R[:,0],self.fsam_R[:,2],color=self.pp.Rcolor,alpha=0.2)
        ax.axhline(self.Rmean,ls="--",lw=2,color=self.pp.Rcolor)
        ax.plot(self.phsams_nq,self.fsam_nq[:,1],color=self.pp.Qcolor,label="Moisture convergence (Q)")
        ax.fill_between(self.phsams_nq,self.fsam_nq[:,0],self.fsam_nq[:,2],color=self.pp.Qcolor,alpha=0.2)
        ax.axhline(self.Qmean,ls="--",lw=2,color=self.pp.Qcolor)
        ax.legend(loc=(0.0,1),ncol=4,frameon=False,fontsize=10)    
        ax.set_ylabel(f"""Exchange
    ({self.pp.uflux})""",fontsize=12)
        ax.margins(0)
        ymin,ymax=ax.get_ylim()
        ax.set_ylim(ymin,ymax*1.2)
        ax.set_xlabel(f"""Phase""",fontsize=12)
        ax.text(0.5,0.95,f"""{self.name}""",va="top",ha="center",fontsize=12,transform=ax.transAxes)
        ax.text(0.5,0.90,rf"""Global imbalance $\langle Q\rangle-\langle R\rangle$ = {self.Qmean-self.Rmean:+.3f}$\times 10^3$ km$^3$/yr""",va="top",ha="center",fontsize=8,transform=ax.transAxes)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.0,hspace=0.15)    