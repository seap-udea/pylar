#############################################################
# REQUIRED PACAKGES
#############################################################
from ipylar import *
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from dateutil.relativedelta import relativedelta
import matplotlib.dates as mdates

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
        vnq=0.05 #Molina, according to analysis of uncertainties in data processing of ERA5
        data_atm["dnq"]=vnq*abs(data_atm["nq_mass"])/1e3

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
        vR=0.15 
        data_riv["dR"]=vR*abs(data_riv["R"])

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
        self.tini=float(data_merge["time"].iloc[0])
        self.tend=float(data_merge["time"].iloc[-1])
        self.dini=data_merge["datetime"].iloc[0]
        self.tspan=data_merge["time"].iloc[-1]-data_merge["time"].iloc[0]
        self.yspan=self.tspan/Const.years
        self.dspan=self.tspan/Const.days

        #Global properties
        self.ndata=len(data_merge)
        self.nmerge=len(data_merge)
        self.vnq=vnq
        self.vR=vR

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

        #Interpolation
        self.Rint=interp1d(self.data["time"],self.data["R"],kind='slinear')
        self.nqint=interp1d(self.data["time"],self.data["nq"],kind='slinear')
            
        #Data global properties
        print(f"*"*40)
        print(f"{self.key}")
        print(f"Number of data points: merge = {self.nmerge}, total = {self.ndata}")
        print(f"Initial time: {self.tini}")
        print(f"Initial date: {self.dini}")
        print(f"Total data points: {self.ndata}")
        print(f"Time span: {self.yspan} a = {self.dspan} d")
        print(f"Signal correlation (R - nq): {self.rhoRnq}")

    def plot_basin_series(
            self,
            ti_win=2000,dt_win=5,
            over_interp=10,
            show_ref=False):

        # Properties
        #     
        #Safe colors
        Qcolor="#d95f02"
        Rcolor="#7570b3"

        #Common strings
        nam_nablaq="Moisture convergence"
        str_nablaq=r"$Q$"
        nam_R="River discharge"
        str_R=r"$R$"
        sty_nq=dict(color=Rcolor,alpha=1)
        sty_R=dict(color=Qcolor,alpha=1)
        lab_nq=rf"{nam_nablaq} ({str_nablaq})"
        lab_R=rf"{nam_R} ({str_R})"

        #Units
        nam_flu="Flux"
        UQ=1e9/(365.25*86400) #dam^3 -> m^3
        uflux=r"$10^3$ km$^3$/yr"
        #Time
        utime="d"
        UT=1
        #Total stored
        UI=1e3 #km^3 -> 10^3 km^3

        #Full signal
        tfulls=np.linspace(self.data["time"].iloc[0],self.data["time"].iloc[-1],over_interp*self.ndata)
        dfulls=[datetime.fromtimestamp(t) for t in tfulls]
        Rfulls=self.Rint(tfulls)
        nqfulls=self.nqint(tfulls)

        #Plot full signal
        fig,axs=plt.subplots(2,1,figsize=(9,8))

        #Full panel
        ax=axs[0]

        ax.plot(self.data["datetime"],self.data["nq"]/UQ/UI,label=lab_nq,**sty_nq)
        ax.plot(self.data["datetime"],self.data["R"]/UQ/UI,label=lab_R,**sty_R)

        if show_ref:
            ax.plot(self.data["datetime"],self.data["nqref"]/UQ/UI,label=f"{str_nablaq} ref.",**sty_nq,ls=':')
            ax.plot(self.data["datetime"],self.data["R"]/UQ/UI,label=f"{str_R} ref.",**sty_R,ls=":")
            
        ax.set_ylabel(f"Flux [{uflux}]");
        legend=ax.legend(loc=(0.0,1),ncol=6,frameon=False);

        #Zoom
        ax=axs[1]
        ax.plot(dfulls,nqfulls/UQ/UI,label=rf"{nam_R} ({str_R})",**sty_nq)
        ax.errorbar(self.data["datetime"],self.data["nq"]/UQ/UI,yerr=self.data["dnq"]/UQ/UI,ls='None',capsize=5,**sty_nq)
        ax.plot(dfulls,Rfulls/UQ/UI,label=rf"{nam_R} ({str_R})",**sty_R)
        ax.errorbar(self.data["datetime"],self.data["R"]/UQ/UI,yerr=self.data["dR"]/UQ/UI,ls='None',capsize=5,**sty_R)
        
        if show_ref:
            #Interpolation
            self.Rint=interp1d(self.data["time"],self.data["Rref"],kind='slinear')
            self.nqint=interp1d(self.data["time"],self.data["nqref"],kind='slinear')

            #Full signal
            tfulls=np.linspace(self.data["time"].iloc[0],self.data["time"].iloc[-1],over_interp*self.ndata)
            dfulls=[datetime.fromtimestamp(t) for t in tfulls]
            Rfulls=self.Rint(tfulls)
            nqfulls=self.nqint(tfulls)

            ax.plot(dfulls,nqfulls/UQ/UI,**sty_nq,ls=":")
            ax.errorbar(self.data["datetime"],self.data["nqref"]/UQ/UI,yerr=self.data["dnq"]/UQ/UI,ls='None',capsize=5,**sty_nq)
            ax.plot(dfulls,Rfulls/UQ/UI,**sty_R,ls=":")
            ax.errorbar(self.data["datetime"],self.data["Rref"]/UQ/UI,yerr=self.data["dR"]/UQ/UI,ls='None',capsize=5,**sty_R)

        ax.set_ylabel(f"Flux [{uflux}]");
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
