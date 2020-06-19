{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import xarray as xr\n",
    "import numpy as np\n",
    "from matplotlib import pyplot as plt\n",
    "from cartopy import crs as ccrs\n",
    "import cartopy.feature as cfeature\n",
    "from cartopy.util import add_cyclic_point"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "def graph_year(ds, date, plot, title, coords):\n",
    "        \n",
    "    clevs = np.array([0,1,2,3,4,120])\n",
    "    colorange = ['red', 'orange', 'yellow','green','blue']\n",
    "    crs = ccrs.PlateCarree()\n",
    "    X = ds['xt_ocean']\n",
    "    Y = ds['yt_ocean']\n",
    "    if type(date) == str:\n",
    "        Z = ds['MI'].sel(time=date).squeeze()\n",
    "    if type(date) == int:\n",
    "        Z = ds['MI'].sel(year=date).squeeze()\n",
    "    Z, X = add_cyclic_point(Z,coord=X)\n",
    "    im = plot.contourf(X,Y,Z,clevs,colors=colorange,transform=crs)\n",
    "    if plot == axs[0,0]:\n",
    "        cbar = plt.colorbar(im,ax=axs,orientation='horizontal',fraction=0.05,pad=0.05,shrink=0.5)\n",
    "        cbar.set_label('metabolic index',fontsize=12)\n",
    "    \n",
    "    # Zoom in on a region\n",
    "    if coords != None:\n",
    "        ax.set_extent(coords)\n",
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
