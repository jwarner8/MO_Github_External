"""
Script to plot the 4p4km NWP fields to identify drylines, supporting the IOP plane obs
"""

__author__ = "James Warner"
__email__  = "james.warner@metoffice.gov.uk"

from matplotlib import pyplot as plt
import matplotlib
import matplotlib.ticker as mticker
matplotlib.use('Agg')
import cartopy.crs as ccrs
import cartopy
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import os
import subprocess
import iris
import datetime as dt

# Location to store data and plots
data_loc = os.environ['SCRATCH']+'/0201_kapex/'


def get_data(init):

    """
    Get required fields from MASS, and store on $SCRATCH
    """

    print('Creating filtered retrieval file for',init)
    queryfile = open(data_loc+'/qf', 'w')
    queryfile.write('begin\n')
    queryfile.write('stash=(3225,3226,3236,3250,2205,16222,3245,4201,16202,16204)'+'\n')
    queryfile.write('pp_file=(')

    for leadtime in range(0,37,12):

        queryfile.write('"'+init+'_kapex_4p4km_RAL3p2_pvera'+str(leadtime).zfill(3)+'.pp", ')
        queryfile.write('"'+init+'_kapex_4p4km_RAL3p2_pverb'+str(leadtime).zfill(3)+'.pp", ')
        queryfile.write('"'+init+'_kapex_4p4km_RAL3p2_pverd'+str(leadtime).zfill(3)+'.pp", ')

    # Write latest contents to file
    queryfile.close()

    # Remove last ', ' and replace with a bracket/new line to close the queryfile
    with open(data_loc+'/qf', 'rb+') as filehandle:
        filehandle.seek(-2, os.SEEK_END)
        filehandle.truncate()

    queryfile = open(data_loc+'/qf', 'a')
    queryfile.write(')\n')
    queryfile.write('end')
    queryfile.close()

    # Get data from MASS
    print('Pulling data from MASS')
    subprocess.check_output(['/opt/moose-client-wrapper/bin/moo', 'select',data_loc+'/qf',
                          'moose:devfc/u-db313/field.pp/',data_loc])


def plot(init):
  
    """
    Plot retrieved fields
    """

    # Load all cubes
    cubes = iris.load(data_loc+'/'+init+'*.pp')
    print(cubes)


    # iterate over all timeslices, every 3 hours
    for timeslice in range(0,48,3):
        print(timeslice)

        # Plot figure
        gs = gridspec.GridSpec(300,210)
        crs = ccrs.PlateCarree()

        ax = plt.subplot(gs[5:95,5:95],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
        gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
        gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = False ; gl.left_labels = True ; gl.right_labels = False

        plotslice = cubes.extract('air_pressure_at_sea_level')[0][timeslice,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        cs = ax.contourf(plotslice.coord('longitude').points,plotslice.coord('latitude').points,plotslice.data/100,np.linspace(1005,1020,16),transform=crs,extend='both')
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
        ax.set_title('MSLP (hPa)')
        ax.set_extent([9,31,-33,-11])
        ax = plt.subplot(gs[5:95,96:99])
        plt.colorbar(cs,cax=ax, orientation='vertical')

        ax = plt.subplot(gs[5:95,105+10:195+10],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
        gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
        gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = False ; gl.left_labels = False ; gl.right_labels = False

        plotslice = cubes.extract('air_temperature')[0][timeslice,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        cs = ax.contourf(plotslice.coord('longitude').points,plotslice.coord('latitude').points,plotslice.data,np.linspace(285,315,31),transform=crs,extend='both')
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
        ax.set_title('1.5m AirT (K)')
        ax.set_extent([9,31,-33,-11])
        ax = plt.subplot(gs[5:95,197+10:199+10])
        plt.colorbar(cs,cax=ax, orientation='vertical')

        ax = plt.subplot(gs[105:195,5:95],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
        gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
        gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = False ; gl.left_labels = True ; gl.right_labels = False

        plotslice = cubes.extract('dew_point_temperature')[0][timeslice,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        cs = ax.contourf(plotslice.coord('longitude').points,plotslice.coord('latitude').points,plotslice.data,np.linspace(265,295,31),transform=crs,extend='both')
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
        ax.set_title('1.5m DP (K)')
        ax.set_extent([9,31,-33,-11])
        ax = plt.subplot(gs[105:195,96:99])
        plt.colorbar(cs,cax=ax, orientation='vertical')

        ax = plt.subplot(gs[105:195,105+10:195+10],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
        gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
        gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = False ; gl.left_labels = False ; gl.right_labels = False

        plotslice1 = cubes.extract('x_wind')[0][timeslice,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        plotslice2 = cubes.extract('y_wind')[0][timeslice,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        plotslicemag = np.sqrt((plotslice1.data**2)+(plotslice2.data**2))
        cs = ax.contourf(plotslice1.coord('longitude').points,plotslice1.coord('latitude').points,plotslicemag,np.linspace(0,10,21),cmap='Reds',transform=crs,extend='both')
        ax.streamplot(plotslice1.coord('longitude').points,plotslice1.coord('latitude').points,plotslice1.data,plotslice2.data,linewidth=0.4,density=1.5,color='black',transform=ccrs.PlateCarree())
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
        ax.set_title('10m Wind (m/s)')
        ax.set_extent([9,31,-33,-11])
        ax = plt.subplot(gs[105:195,196+10:199+10])
        plt.colorbar(cs,cax=ax, orientation='vertical')

        ax = plt.subplot(gs[205:295,5:95],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
        gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
        gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = True ; gl.left_labels = True ; gl.right_labels = False

        olrslice = cubes.extract('toa_outgoing_longwave_flux')[0][timeslice,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        cs = ax.contourf(olrslice.coord('longitude').points,olrslice.coord('latitude').points,(olrslice.data/(5.67*10**-8))**0.25,np.linspace(190,290,21),transform=crs,extend='both')
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
        ax.set_extent([9,31,-33,-11])  
        ax.set_title('TOA BT (K)') 
        ax = plt.subplot(gs[205:295,96:99])
        plt.colorbar(cs,cax=ax, orientation='vertical')

        ax = plt.subplot(gs[205:295,105+10:195+10],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
        gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
        gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = True ; gl.left_labels = False ; gl.right_labels = False

        plotslice = cubes.extract('geopotential_height')[0][int(timeslice/3),:,:,:].intersection(longitude=(9,31),latitude=(-33,-11))
        z_700 = plotslice.extract(iris.Constraint(pressure=700))
        z_850 = plotslice.extract(iris.Constraint(pressure=850))
        cs = ax.contourf(plotslice.coord('longitude').points,plotslice.coord('latitude').points,z_700.data-z_850.data,np.linspace(1620,1700,21),cmap='terrain',transform=crs,extend='both')
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
        ax.set_extent([9,31,-33,-11])
        ax.set_title('z700-z850 (m)')
        ax = plt.subplot(gs[205:295,196+10:199+10])
        plt.colorbar(cs,cax=ax, orientation='vertical')

        inittime  = dt.datetime.strptime(init,'%Y%m%dT%H%MZ')
        validtime = dt.datetime.strptime(init,'%Y%m%dT%H%MZ') + dt.timedelta(hours=timeslice)

        plt.suptitle('MetUM KAPEX 4p4km  |  Init '+inittime.strftime('%d/%m %H%MZ')+' | Valid '+validtime.strftime('%d/%m %H%MZ'),weight='bold',y=0.91)
        fig = plt.gcf()
        fig.set_size_inches(8,11)
        plt.savefig(data_loc+'/'+init+'_'+str(timeslice).zfill(2)+'.png',dpi=150,bbox_inches='tight')
        plt.close()


def plot_dryline_structure():

    """
    For flight planning, look at gradient of relative humidity at different pressure
    levels
    """

    for timeslice in range(0,17):

        cube = iris.load_cube(data_loc+'/20240112T0000Z*verd*.pp','relative_humidity')[timeslice,-4,:,:]
        print(cube)

        inittime  = dt.datetime.strptime('20240112T0000Z','%Y%m%dT%H%MZ')
        validtime = dt.datetime.strptime('20240112T0000Z','%Y%m%dT%H%MZ') + dt.timedelta(hours=timeslice*3)

        gs = gridspec.GridSpec(100,100)
        crs = ccrs.PlateCarree()

        ax = plt.subplot(gs[5:95,5:95],projection=crs)
        ax.add_feature(cartopy.feature.COASTLINE, linestyle='-')
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
        gl = ax.gridlines(draw_labels=True,linestyle='--')
    #    gl.xlocator = mticker.FixedLocator([10,14,18,22,26,30])
    #    gl.ylocator = mticker.FixedLocator([-32,-28,-24,-20,-16,-12])
        gl.top_labels = False ; gl.bottom_labels = True ; gl.left_labels = True ; gl.right_labels = False

        cs = ax.contourf(cube.coord('longitude').points,cube.coord('latitude').points,cube.data,np.linspace(0,100,21),cmap='terrain',transform=crs,extend='both')
        plt.plot(20.6,-26.53,'x',markersize=6,color='red',transform=crs)
        plt.plot(22.61,-26.34,'x',markersize=6,color='red',transform=crs)
        plt.plot(20,-26.44,'x',markersize=6,color='red',transform=crs)
      #  plt.plot(21.3,-28.4,'o',markersize=8,color='black',transform=crs)
        ax.set_title('MetUM KAPEX 4p4km  |  Init '+inittime.strftime('%d/%m %H%MZ')+' | Valid '+validtime.strftime('%d/%m %H%MZ'))
        ax.set_extent([16,26,-31,-21])
        ax = plt.subplot(gs[5:95,96:99])
        plt.colorbar(cs,cax=ax, label='RelH 850hPa (%)',orientation='vertical')

        fig = plt.gcf()
        fig.set_size_inches(6,6)
        plt.savefig('temp_'+str(timeslice).zfill(2)+'.png',dpi=150,bbox_inches='tight')
        plt.close()


def main():

#    # Set model forecast
    init = '20240115T0000Z'

    # Get data
    get_data(init)

    # Plot data
    plot(init)
    subprocess.check_output(['/usr/bin/convert', '-loop', '0', '-delay', '180', data_loc+'/'+init+'*.png', data_loc+'/'+init+'.gif'])
    subprocess.check_output(['/usr/bin/ffmpeg', '-y','-i', data_loc+'/'+init+'.gif', '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', data_loc+'/'+init+'.mp4'])

    # Plot relative humidity spatially.
#    plot_dryline_structure()
#    subprocess.check_output(['/usr/bin/convert', '-loop', '0', '-delay', '180', 'temp*.png', 'temp.gif'])
#    subprocess.check_output(['/usr/bin/ffmpeg', '-y','-i', 'temp.gif', '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', 'temp.mp4'])  


if __name__ == '__main__':

    main()
