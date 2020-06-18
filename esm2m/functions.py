{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "import xarray as xr\n",
    "from matplotlib import pyplot as plt\n",
    "from cartopy import crs as ccrs\n",
    "import cartopy.feature as cfeature\n",
    "from cartopy.util import add_cyclic_point\n",
    "import numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Load the area data\n",
    "rootdir = '/local/ss23/GFDL_LEs/'\n",
    "subdir = 'AREA_FILES_ETC'\n",
    "filename_area = 'WOA2001_grid.nc'\n",
    "path_area = rootdir+subdir+'/'+filename_area\n",
    "area = xr.open_dataset(path_area)['AREA'].rename({'latitude':'yt_ocean','longitude':'xt_ocean'})"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [],
   "source": [
    "def graph(ds, date, plot, title, coords):\n",
    "        \n",
    "    clevs = np.array([0,1,2,3,4,120])\n",
    "    colorange = ['red', 'orange', 'yellow','green','blue']\n",
    "    crs = ccrs.PlateCarree()\n",
    "    X = ds['xt_ocean']\n",
    "    Y = ds['yt_ocean']\n",
    "    \n",
    "    if type(date) == str:\n",
    "        Z = ds['MI'].sel(time=date).squeeze()\n",
    "    if type(date) == int:\n",
    "        Z = ds['MI'].sel(year=date).squeeze()\n",
    "    Z, X = add_cyclic_point(Z,coord=X)\n",
    "    im = plot.contourf(X,Y,Z,clevs,colors=colorange,transform=crs)\n",
    "\n",
    "    if coords != None:\n",
    "        plot.set_extent(coords)\n",
    "\n",
    "    # Add a land mask to your plot, as well as grid lines and coastlines\n",
    "    plot.add_feature(cfeature.LAND,zorder=10,facecolor='darkgray')\n",
    "    plot.gridlines()\n",
    "    plot.coastlines()\n",
    "    plot.set_title(title,fontsize=14,loc='center')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "def year_comp(ds, plot, minLat, maxLat, minLon, maxLon, cont1, cont2, cont3):\n",
    "    area = area.assign_coords({'xt_ocean':ds['xt_ocean'],'yt_ocean':ds['yt_ocean']})\n",
    "    years = np.arange(1990, 2110, 10)\n",
    "    colors = ['khaki','limegreen','orange','teal','mediumpurple','slategray','indigo','palegreen','turquoise','pink','navy']\n",
    "    colorblue = ['turquoise','darkcyan','powderblue','deepskyblue','royalblue','navy','slategray','blue']\n",
    "    colortab = ['tab:blue','tab:orange','tab:green','tab:purple','tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','pink','navy','darkgreen']\n",
    "    \n",
    "    if minLat != None:\n",
    "        ds_lim = ds.sel(yt_ocean=slice(minLat,maxLat))\n",
    "    if minLon != None:\n",
    "        ds_lim = ds.sel(yt_ocean=slice(minLat,maxLat))\n",
    "    oceanmask = np.isfinite(ds_z0avg['MI'].sel(time='1950-01-16',yt_ocean=slice(minLat,maxLat)).squeeze())\n",
    "    area_lim = area.sel(yt_ocean=slice(minLat,maxLat))\n",
    "    area_masked = area_lim.where(oceanmask,np.nan)\n",
    "    \n",
    "    i = 0\n",
    "    for year in years:\n",
    "        ds_year = ds_lim.sel(time=slice(str(year)+'-01-01',str(year)+'-12-31'))\n",
    "        mean_year = (ds_year['MI']*area_masked).sum(['xt_ocean','yt_ocean'])/(area_masked.sum(['xt_ocean','yt_ocean']))\n",
    "        plot.plot(np.unique(mean_year['time.month']),mean_year.groupby('time.month').mean(),color=colortab[i],label=str(year))\n",
    "        i += 1\n",
    "    plot.legend()\n",
    "    if cont1 == True:\n",
    "        plot.axhline(y=1.0, xmin=0,xmax=1,color='red',dashes=[6,2])\n",
    "    if cont2 == True:\n",
    "        plot.axhline(y=2.0, xmin=0,xmax=1,color='orange',dashes=[6,2])\n",
    "    if cont3 == True:\n",
    "        plot.axhline(y=3.0, xmin=0,xmax=1,color='yellow',dashes=[6,2])\n",
    "    plot.set_ylabel('metabolic index')\n",
    "    plot.set_xlabel('month')\n",
    "    plot.set_title('Mean Metabolic Index, '+str(minLat)+ ' to '+str(maxLat)+ ', by Year')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "def graph_ens(ds, title, date, coords):\n",
    "    fig, axs = plt.subplots(figsize=(20,20),nrows=6,ncols=5, subplot_kw={'projection':ccrs.Robinson(central_longitude=180)})\n",
    "    fig.suptitle(title,fontsize=20)\n",
    "    wn.filterwarnings('ignore')\n",
    "    ensNum = 0\n",
    "    for row in range(6):\n",
    "        for col in range(5):\n",
    "            ds_ens = ds.sel(ensemble=ensNum)\n",
    "            graph(ds_ens,date,axs[row,col],'Ensemble #'+str(ensNum+1), coords)\n",
    "            ensNum += 1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.2"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
