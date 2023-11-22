class Plotting:
    def plot_data(self, out_file, data):
        """Plots the seasonal_means to screen."""
        # Plot the data:
        X, Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')

        # Define plot parameters:
        nplots = 3
        max_val = [80, 80, 20]
        nbound = [21, 21, 21]
        unit = ['[pptv]','[pptv]','unitless']
        plot_title = ['TROPOMI cloud-sliced NO2',
                      'TROPOMI cloud-sliced NO2 error',
                      'Number of cloud-sliced data points']

        fig, ax = plt.subplots(3, 1, figsize=(5,11), subplot_kw=dict(projection=ccrs.PlateCarree()))

        # Plot the subplots:
        for i in range(nplots):

            # Define what data to plot:
            if i==0: plot_vals = np.squeeze(self.mean_gno2vmr)
            if i==1: plot_vals = np.squeeze(self.mean_gerr)
            if i==2: plot_vals = np.squeeze(self.gcnt)

            if i==0: 
                tickval = [0,25,50,75,100]
            else:
                tickval = [0,5,10,15,20]

            ax[i].coastlines(resolution='50m')

            ax[i].set_extent([-179.9, 179.9, -75, 75], crs=ccrs.PlateCarree())

            data_crs = ccrs.PlateCarree()
            bounds = np.linspace(0, max_val[i],nbound[i])

            c = ax[i].pcolormesh(X,Y, plot_vals, transform=data_crs,cmap='jet', 
                                 vmin=0, vmax=max_val[i])

            cb = fig.colorbar(c, ax=ax[i],label=unit[i], orientation='horizontal',
                              shrink=0.5,pad=0.01,boundaries=bounds,ticks=tickval )
            
            cb.ax.tick_params(labelsize=10, direction='in', length=6)

        plt.savefig(out_file, format='png')
        plt.show()