import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from cartopy import crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import warnings as wn

def graph(ds, plot, title, date, coords, ens):
        
    clevs = np.array([0,1,2,3,4,120])
    colorange = ['red', 'orange', 'yellow','green','blue']
    crs = ccrs.PlateCarree()
    X = ds['xt_ocean']
    Y = ds['yt_ocean']
    if date == None:
        Z = ds['MI'].squeeze()
    elif type(date) == str:
        Z = ds['MI'].sel(time=date).squeeze()
    elif type(date) == int:
        Z = ds['MI'].groupby('time.year').mean().sel(year=date).squeeze()
    Z, X = add_cyclic_point(Z,coord=X)
    im = plot.contourf(X,Y,Z,clevs,colors=colorange,transform=crs)
    
    # Zoom in on a region
    if coords != None:
        plot.set_extent(coords)

    # Add a land mask to your plot, as well as grid lines and coastlines
    plot.add_feature(cfeature.LAND,zorder=10,facecolor='darkgray')
    plot.gridlines()
    plot.coastlines()
    plot.set_title(title,fontsize=14,loc='center')
    
    if ens:
        return im
    cbar = plt.colorbar(im,ax=plot,orientation='horizontal',fraction=0.05,pad=0.05,shrink=0.5)
    cbar.set_label('metabolic index',fontsize=12)
    
    
def graph_ens(ds, title, date, coords):
    fig, axs = plt.subplots(figsize=(20,20),nrows=6,ncols=5, subplot_kw={'projection':ccrs.Robinson(central_longitude=180)})
    fig.suptitle(title,fontsize=20)
    wn.filterwarnings('ignore')
    ensNum = 0
    for row in range(6):
        for col in range(5):
            ds_ens = ds.sel(ensemble=ensNum)
            im = graph(ds_ens, axs[row,col], 'Ens '+str(ensNum+1), date, coords, True)
            ensNum += 1
    cbar = plt.colorbar(im,ax=axs,orientation='horizontal',fraction=0.05,pad=0.05,shrink=0.5)
    cbar.set_label('metabolic index',fontsize=12)
    
def percent_plot(ds, plot, title, red, orange, yellow, col, lb, depth, total):
    if red:
        ds_red = ~np.isnan(ds.where(ds['MI']<1.0))
        ds_rSum = ds_red['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
        ds_rPercent = (ds_rSum/total)*100
        if col == None:
            # plot.plot(np.unique(ds_rPercent['time.year']),ds_rPercent.groupby('time.year').mean(),color='teal',label='MI < 1')
            plot.plot(ds_rPercent['time'],ds_rPercent,color='teal',label='MI < 1')
        else:
            # plot.plot(np.unique(ds_rPercent['time.year']),ds_rPercent.groupby('time.year').mean(),color=col,label=lb)
            plot.plot(ds_rPercent['time'],ds_rPercent,color=col,label=lb)
        if depth == 'k01':
            plot.set_ylim(-0.2,6)
        if depth == 'irr0.1' or depth == 'k11':
            plot.set_ylim(0,1.5)
        
    if orange:
        ds_orange = ~np.isnan(ds.where(ds['MI']<2.0))
        ds_oSum = ds_orange['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
        ds_oPercent = (ds_oSum/total)*100
        if col == None:
            # plot.plot(np.unique(ds_oPercent['time.year']),ds_oPercent.groupby('time.year').mean(),color='dodgerblue',label='MI < 2')
            plot.plot(ds_oPercent['time'],ds_oPercent,color='dodgerblue',label='MI < 2')
        else:
            # plot.plot(np.unique(ds_oPercent['time.year']),ds_oPercent.groupby('time.year').mean(),color=col,label=lb)
            plot.plot(ds_oPercent['time'],ds_oPercent,color=col,label=lb)
        if depth == 'k01':
            plot.set_ylim(15,30)
        if depth == 'irr0.1':
            plot.set_ylim(3,11)
        if depth == 'k11':
            plot.set_ylim(2,15)
    
    if yellow:
        if col == None:
            col = 'cyan'
        if lb == None:
            lb = 'MI < 3'
        ds_yellow = ~np.isnan(ds.where(ds['MI']<3.0))
        ds_ySum = ds_yellow['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
        ds_yPercent = (ds_ySum/total)*100
        if col == None:
            # plot.plot(np.unique(ds_yPercent['time.year']),ds_yPercent.groupby('time.year').mean(),color='cyan', label='MI < 3')
            plot.plot(ds_yPercent['time'],ds_yPercent,color='cyan',label='MI < 3')
        else:
            # plot.plot(np.unique(ds_yPercent['time.year']),ds_yPercent.groupby('time.year').mean(),color=col, label=lb)  
            plot.plot(ds_yPercent['time'],ds_yPercent,color=col,label=lb)
        if depth == 'k01':
            plot.set_ylim(28,40)
        if depth == 'irr0.1':
            plot.set_ylim(12,28)
        if depth == 'k11':
            plot.set_ylim(15,27)
        
    if red and yellow:
        if depth == 'k01':
            plot.set_ylim(-1,40)
        if depth == 'irr0.1':
            plot.set_ylim(0,24)
        if depth == 'k11':
            plot.set_ylim(0,28)
        plot.legend()
    elif orange and yellow:
        plot.set_ylim(14,40)
        plot.legend()
    elif red and orange:
        plot.set_ylim(-1,30)
        plot.legend()
        
    plot.set_title(title)
    plot.set_ylabel('Percent of Ocean')
    
def percent_ens(ds, title, red, yellow, orange, depth, total):
    ds_lim = ds.sel(time=slice('1990-01-16','2100-12-16'))
    fig, axs = plt.subplots(figsize=(20,20),nrows=6,ncols=5)
    fig.suptitle(title,fontsize=20) # Specify a figure title
    wn.filterwarnings('ignore')
    ensNum = 0
    for row in range(6):
        for col in range(5): 
            percent_plot(ds_lim.sel(ensemble=ensNum), axs[row,col], 'Ens '+str(ensNum+1), red, yellow, orange, None, None, depth, total)
            ensNum += 1

            
def year_comp(area, ds, plot, coords, cont1, cont2, cont3, depth):
    
    years = np.arange(1990, 2110, 10)
    colors = ['khaki','limegreen','orange','teal','mediumpurple','slategray','indigo','palegreen','turquoise','pink','navy']
    colorblue = ['turquoise','darkcyan','powderblue','deepskyblue','royalblue','navy','slategray','blue']
    colortab = ['tab:blue','tab:orange','tab:green','tab:purple','tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','pink','navy','darkgreen']
    
    ds_lim = ds.sel(xt_ocean=slice(coords[0],coords[1]),yt_ocean=slice(coords[2],coords[3]))
    oceanmask = np.isfinite(ds['MI'].sel(time='1950-01-16',xt_ocean=slice(coords[0],coords[1]),yt_ocean=slice(coords[2],coords[3])).squeeze())
    area_lim = area.sel(xt_ocean=slice(coords[0],coords[1]),yt_ocean=slice(coords[2],coords[3]))
        
    area_masked = area_lim.where(oceanmask,np.nan)
    i = 0
    for year in years:
        ds_year = ds_lim.sel(time=slice(str(year)+'-01-01',str(year)+'-12-31'))
        mean_year = (ds_year['MI']*area_masked).sum(['xt_ocean','yt_ocean'])/(area_masked.sum(['xt_ocean','yt_ocean']))
        plot.plot(np.unique(mean_year['time.month']),mean_year.groupby('time.month').mean(),color=colortab[i],label=str(year))
        i += 1
    plot.legend()
    if cont1 == True:
        plot.axhline(y=1.0, xmin=0,xmax=1,color='red',dashes=[6,2])
    if cont2 == True:
        plot.axhline(y=2.0, xmin=0,xmax=1,color='orange',dashes=[6,2])
    if cont3 == True:
        plot.axhline(y=3.0, xmin=0,xmax=1,color='yellow',dashes=[6,2])
    plot.set_ylabel('metabolic index')
    plot.set_xlabel('month')
    if coords[0] == None:
        plot.set_title('Mean Metabolic Index, latitude '+str(coords[2])+ ' to '+str(coords[3])+ ', by Year, '+depth)
    elif coords[2] == None:
        plot.set_title('Mean Metabolic Index, longitude'+str(coords[0])+ ' to '+str(coords[1])+ ', by Year, '+depth)
    else:
        plot.set_title('Mean Metabolic Index, ['+str(coords[0])+','+str(coords[1])+','+str(coords[2])+','+str(coords[3])+'] by Year, '+depth)
        
def percent_by_year(ds, title, plot, total):
    colortab = ['tab:blue','tab:orange','tab:green','tab:purple','tab:brown','pink','tab:grey','tab:olive','tab:cyan','tab:pink','navy','darkgreen']
    years = np.arange(1990, 2110, 10)
    i = 0
    for year in years:
        ds_year = ds.sel(time=slice(str(year)+'-01-01',str(year)+'-12-31'))
#         ds_red = ~np.isnan(ds_year.where(ds_year['MI']<mi))
#         ds_rSum = ds_red['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
#         ds_rPercent = (ds_rSum/total)*100
        plot.plot(ds_year,color=colortab[i],label=str(year))
        i += 1
    # plot.legend()
    plot.set_title(title)
    plot.set_xlabel('months')
    plot.set_ylabel('percent of ocean')

def percent_ens_comp(ds_all, plot, mi, total):
    colorlight = ['skyblue', 'bisque', 'lightgreen', 'thistle', 'peachpuff', 'lightpink', 'gainsboro', 'palegoldenrod', 'lightcyan', 'lavenderblush', 'lightsteelblue', 'palegreen']
    years = np.arange(1990, 2110, 10)
    ens = np.arange(0,30,1)
    i = 0
    wn.filterwarnings('ignore')
    for mem in ens:
        ds = ds_all.sel(ensemble=mem)
        i = 0
        for year in years:
            ds_year = ds.sel(time=slice(str(year)+'-01-01',str(year)+'-12-31'))
            ds_red = ~np.isnan(ds_year.where(ds_year['MI']<mi))
            ds_rSum = ds_red['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
            ds_rPercent = (ds_rSum/total)*100
            plot.plot(ds_rPercent,color=colorlight[i])
            i += 1
    
def months_of_year(data,year,title):
    fig, axs = plt.subplots(figsize=(20,20),nrows=4,ncols=3, subplot_kw={'projection':ccrs.Robinson(central_longitude=180)})
    fig.suptitle(title,fontsize=15) # Specify a figure title
    graph(data.sel(time=slice(year+'-01-01',year+'-01-31')), axs[0,0], '(a) January', None, None, True)
    graph(data.sel(time=slice(year+'-02-01',year+'-02-28')), axs[0,1], '(b) February', None, None, True)
    graph(data.sel(time=slice(year+'-03-01',year+'-03-31')), axs[0,2], '(c) March', None, None, True)
    graph(data.sel(time=slice(year+'-04-01',year+'-04-30')), axs[1,0], '(d) April', None, None, True)
    graph(data.sel(time=slice(year+'-05-01',year+'-05-31')), axs[1,1], '(e) May', None, None, True)
    graph(data.sel(time=slice(year+'-06-01',year+'-06-30')), axs[1,2], '(f) June', None, None, True)
    graph(data.sel(time=slice(year+'-07-01',year+'-07-31')), axs[2,0], '(g) July', None, None, True)
    graph(data.sel(time=slice(year+'-08-01',year+'-08-31')), axs[2,1], '(h) August', None, None, True)
    graph(data.sel(time=slice(year+'-09-01',year+'-09-30')), axs[2,2], '(e) September', None, None, True)
    graph(data.sel(time=slice(year+'-10-01',year+'-10-30')), axs[3,0], '(f) October', None, None, True)
    graph(data.sel(time=slice(year+'-11-01',year+'-11-30')), axs[3,1], '(g) November', None, None, True)
    im = graph(data.sel(time=slice(year+'-12-01',year+'-12-31')), axs[3,2], '(h) December', None, None, True)
    cbar = plt.colorbar(im,ax=axs,orientation='horizontal',fraction=0.05,pad=0.05,shrink=0.5)
    cbar.set_label('metabolic index',fontsize=12)

def get_red(ds, total):
    ds_red = ~np.isnan(ds.where(ds['MI']<1.0))
    ds_rSum = ds_red['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
    ds_rPercent = (ds_rSum/total)*100
    ds_rAvg = ds_rPercent.mean(dim='ensemble')
    return ds_rAvg

def get_orange(ds, total):
    ds_orange = ~np.isnan(ds.where(ds['MI']<2.0))
    ds_oSum = ds_orange['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
    ds_oPercent = (ds_oSum/total)*100
    ds_oAvg = ds_oPercent.mean(dim='ensemble')
    return ds_oAvg

def get_yellow(ds, total):
    ds_yellow = ~np.isnan(ds.where(ds['MI']<3.0))
    ds_ySum = ds_yellow['MI'].sum(dim='xt_ocean').sum(dim='yt_ocean')
    ds_yPercent = (ds_ySum/total)*100
    ds_yAvg = ds_yPercent.mean(dim='ensemble')
    return ds_yAvg