"""
Calculate apparent heat source, cited from Ye et al. 2023.
"""

__author__ = "James Warner"
__email__  = "james.warner@metoffice.gov.uk"

import subprocess
import os
import glob
import datetime as dt
import iris
iris.FUTURE.datum_support = True
import iris.cube
from iris.coords import DimCoord
import numpy as np
import sys


def calc_appheatsource(t,u,v,w,z,strid,omega_flag=False):

    """
    Calculate apparent heat source given 4D cubes of temperature (t) and 
    winds (u,v,w), with minimum 2 timesteps (needs to calculate a time
    tendency). All pressure levels are assumed to increase from index 0
    (i.e. 200,250,300...).

    Takes approximately 1 hour for 1 degree, with 2 timesteps. 2GB is
    sufficient.

    omega is a bool as to whether w is in m/s (false) or pa/s (true).
    """
 
    print('Running AHS')

    # Check that all cubes have same dimension shape, if not, quit.
    if not all(x == t.shape for x in [t.shape,u.shape,v.shape,w.shape,z.shape]):
        print('Different dims... exiting')
        quit()

    # Get number of time points, levels, lats, lons
    ntime = len(t.coord('time').points)
    nlevs = len(t.coord('pressure').points)
    nlats = len(t.coord('latitude').points)
    nlons = len(t.coord('longitude').points)
    print('Shape ',ntime,nlevs,nlats,nlons)

    # Get pressure levels in Pa.
    plevs = t.coord('pressure').points*100.
    
    # Determine orientation of pressure coord (so appropriate dtheta/dp).
    # Increasing means index 0 200, 1 is 250, 2 is 300 etc.
    if plevs[0] > plevs[1]:
        plev_dir = 'decreasing, not supported'
        quit()

    # Generate cube to save output to nc
    print('Creating output cube')
    diag_cube = iris.cube.Cube(np.zeros((4,ntime-1,nlevs-1,nlats-1,nlons-1)),long_name='absolute_heat_source',
                               dim_coords_and_dims=[(DimCoord([0,1,2,4],long_name='ahs_term'),0),(t[-1].coord('time'),1),
                               (t[:,:-1,:,:].coord('pressure'),2),(t[:,:,1:,:].coord('latitude'),3),(t[:,:,:,1:].coord('longitude'),4)])

    # Go through each gridpoint sequentially and calculate
    for lat in range(1,nlats):
        print(lat,'/',nlats,'|',dt.datetime.now())
        for lon in range(1,nlons):
            for time in range(1,ntime):
                for lev in range(0,nlevs-1):

                    # Calculate local dT/dt, assuming 3 hourly data
                    dTdt = (t[time,lev,lat,lon].data - t[time-1,lev,lat,lon].data)/(3*60*60)

                    # Assuming grid resolution of 1 degree, and approximating to be 111 km, in metres so units match.
                    VdotdelT = u[time,lev,lat,lon].data*((t[time,lev,lat,lon].data-t[time,lev,lat,lon-1].data)/(111*1000))+\
                               v[time,lev,lat,lon].data*((t[time,lev,lat,lon].data-t[time,lev,lat-1,lon].data)/(111*1000))

                    # Calculate potential temperature, and also level below (here lev+1 as increasing/starting at toa).
                    potT = t[time,lev,lat,lon].data * ((100000/plevs[lev])**(287/1004))
                    potT_levbelow = t[time,lev+1,lat,lon].data * ((100000/plevs[lev+1])**(287/1004))

                    if omega_flag == False:
                        # Calculate omega (pa/s) from vertical velocity and geopotential height. Should be -ve when m/s +ve
                        omega = w[time,lev,lat,lon].data*((plevs[lev]-plevs[lev+1])/(z[time,lev,lat,lon].data - z[time,lev+1,lat,lon].data))
                    else:
                        omega = w[time,lev,lat,lon].data

                    # Calculate final term (pressure ratio times w times dTheta/dp (in Pa).
                    ppsdotpotTdw = ((plevs[lev]/100000)**(287/1004))*omega*((potT-potT_levbelow)/(plevs[lev]-plevs[lev+1]))

                    # Total apparent heat source, and components
                    diag_cube.data[0,time-1,lev,lat-1,lon-1] = dTdt
                    diag_cube.data[1,time-1,lev,lat-1,lon-1] = VdotdelT
                    diag_cube.data[2,time-1,lev,lat-1,lon-1] = ppsdotpotTdw
                    diag_cube.data[3,time-1,lev,lat-1,lon-1] = 1004 * (dTdt + VdotdelT + ppsdotpotTdw)
                   
    # Save file.
    iris.save(diag_cube,os.environ['SCRATCH']+'/appheatsource/TMP_ahs_'+strid+'.nc')


def preprocess_um_data(timeid,model):

    """
    An iterator function to break up cubes/pass to absolute heat source
    computation. TimeID is number between 2 and ntime (321), model is 
    str for model configuration. SBATCH requirements in calc_appheatsource
    docstring.
    """

    # Load all cubes for respective model.
    print('Loading cubes')
    all_cubes = iris.load(os.environ['DATADIR']+'/K_SCALE/'+model+'_0p5deg.pp')

    # Get target CTC grid and regrid again to coarse 1 degree.
    if model != "DS_AfricaLAM":
        target_grid = iris.load(os.environ['DATADIR']+'/K_SCALE/DS_RAL3p1_CTC_0p5deg.pp')[0][0,0,:,:].\
                                intersection(longitude=(-30,140),latitude=(-20,30))
    else:
        target_grid = all_cubes[0] # dummy, set as existing LAM domain

    # Regrid to 1 degree resolution using linear regridder.
    lat,lon = target_grid.coord('latitude'),target_grid.coord('longitude')
    lat_min,lon_min = lat.points.min(),lon.points.min()
    lat_max,lon_max = lat.points.max(),lon.points.max()
    latregrd = np.arange(lat_min, lat_max, 1)
    lonregrd = np.arange(lon_min, lon_max, 1)
    target_grid_1deg = target_grid.interpolate([('latitude', latregrd), ('longitude', lonregrd)],
                                                iris.analysis.Linear())
 
    # Extract variables and slice based on timeid
    print('Regridding')
    t = all_cubes.extract('air_temperature')[0][timeid-2:timeid,:,:,:]
    u = all_cubes.extract('x_wind')[0][timeid-2:timeid,:,:,:]
    v = all_cubes.extract('y_wind')[0][timeid-2:timeid,:,:,:]
    w = all_cubes.extract('upward_air_velocity')[0][timeid-2:timeid,:,:,:]
    z = all_cubes.extract('geopotential_height')[0][timeid-2:timeid,:,:,:]

    # Regrid to 1 degree
    t_rgd = t.regrid(target_grid_1deg,iris.analysis.Linear()) 
    u_rgd = u.regrid(target_grid_1deg,iris.analysis.Linear()) 
    v_rgd = v.regrid(target_grid_1deg,iris.analysis.Linear()) 
    w_rgd = w.regrid(target_grid_1deg,iris.analysis.Linear()) 
    z_rgd = z.regrid(target_grid_1deg,iris.analysis.Linear()) 

    # Invoke data so no longer lazy (faster in memory)
    print('Making non-lazy')
    t_rgd.data
    u_rgd.data
    v_rgd.data
    w_rgd.data
    z_rgd.data

    # Call to calculate apparent heat source.
    calc_appheatsource(t_rgd,u_rgd,v_rgd,w_rgd,z_rgd,strid=model+'_T'+str(timeid).zfill(3),omega_flag=False)


def preprocess_era5_data(timeid):

    """
    An iterator function to break up cubes/pass to absolute heat source
    computation. TimeID is number between 2 and ntime (321), model is 
    str for model configuration. SBATCH requirements in calc_appheatsource
    docstring.
    """

    # Get target CTC grid and regrid again to coarse 1 degree.
    target_grid = iris.load(os.environ['DATADIR']+'/K_SCALE/DS_RAL3p1_CTC_0p5deg.pp')[0][0,0,:,:].\
                            intersection(longitude=(-30,140),latitude=(-20,30))

    # Regrid to 1 degree resolution using linear regridder.
    lat,lon = target_grid.coord('latitude'),target_grid.coord('longitude')
    lat_min,lon_min = lat.points.min(),lon.points.min()
    lat_max,lon_max = lat.points.max(),lon.points.max()
    latregrd = np.arange(lat_min, lat_max, 1)
    lonregrd = np.arange(lon_min, lon_max, 1)
    target_grid_1deg = target_grid.interpolate([('latitude', latregrd), ('longitude', lonregrd)],
                                                iris.analysis.Linear())
 
    # Extract variables and slice based on timeid
    print('Loading and Regridding')
    t = iris.load_cube(os.environ['DATADIR']+'/K_SCALE/DS_ERA5_T_0p5deg.nc')[timeid-2:timeid,:,:,:]
    u = iris.load_cube(os.environ['DATADIR']+'/K_SCALE/DS_ERA5_U_0p5deg.nc')[timeid-2:timeid,:,:,:]
    v = iris.load_cube(os.environ['DATADIR']+'/K_SCALE/DS_ERA5_V_0p5deg.nc')[timeid-2:timeid,:,:,:]
    w = iris.load_cube(os.environ['DATADIR']+'/K_SCALE/DS_ERA5_W_0p5deg.nc')[timeid-2:timeid,:,:,:]
    z = iris.load_cube(os.environ['DATADIR']+'/K_SCALE/DS_ERA5_Z_0p5deg.nc')[timeid-2:timeid,:,:,:]

    # Regrid to CTC 1 degree
    t_rgd = t.regrid(target_grid_1deg,iris.analysis.Linear()) 
    u_rgd = u.regrid(target_grid_1deg,iris.analysis.Linear()) 
    v_rgd = v.regrid(target_grid_1deg,iris.analysis.Linear()) 
    w_rgd = w.regrid(target_grid_1deg,iris.analysis.Linear()) 
    z_rgd = z.regrid(target_grid_1deg,iris.analysis.Linear()) 

    # Invoke data so no longer lazy (faster in memory)
    print('Making non-lazy')
    t_rgd.data
    u_rgd.data
    v_rgd.data
    w_rgd.data
    z_rgd.data = z_rgd.data/9.81 # convert geopotential to geopot height.

    # Call to calculate apparent heat source.
    calc_appheatsource(t_rgd,u_rgd,v_rgd,w_rgd,z_rgd,strid='DS_ERA5_T'+str(timeid).zfill(3),omega_flag=True)


def merge_and_integrate(model):

    """
    For some model (passed as string), compute integrated column
    heating, and save as a merged file. Expects all timesteps to
    have been run.
    """

    all_diags = iris.load(os.environ['SCRATCH']+'/appheatsource/TMP_ahs_'+model+'_T*.nc').concatenate()[0]
    print(all_diags)

    # Save all diagnostics prior to column integration, so we can extract profile heating rates.
    print('Saving all merged data')
    iris.save(all_diags, os.environ['SCRATCH']+'/appheatsource/AHS_Merged_fullcolumn_'+model+'.nc')

    # Generate cube to save integrated output to nc
    print('Creating output cube')
    diag_cube = iris.cube.Cube(np.zeros((4,len(all_diags.coord('time').points),len(all_diags.coord('latitude').points),
                               len(all_diags.coord('longitude').points))),long_name='absolute_heat_source',
                               dim_coords_and_dims=[(DimCoord([0,1,2,4],long_name='ahs_term'),0),(all_diags.coord('time'),1),
                               (all_diags.coord('latitude'),2),(all_diags.coord('longitude'),3)])

    # Predefined list of levels - manually defined to include 1000hPa. In SI units.
    levs = [10000,15000,20000,25000,30000,40000,50000,60000,65000,70000,75000,80000,85000,92500,95000,100000]

    # Iterate over all gridpoints and calculate integral.
    for lev in range(0,len(levs)-1):

        # For each plev, integrate Q depending on how large plev layer is (not all equal).
        diag_cube.data[:,:,:,:] = diag_cube.data[:,:,:,:] + (all_diags.data[:,:,lev,:,:]*(levs[lev+1]-levs[lev]))

    # Divide through by gravity
    diag_cube.data[:,:,:,:] = diag_cube.data[:,:,:,:]*(1/9.81)

    # Save output
    iris.save(diag_cube, os.environ['SCRATCH']+'/appheatsource/AHS_Merged_'+model+'.nc')


def main():

    # Preprocess apparent heat source. 6G, 60mins, 2 tasks. 2-321 inc. for full
    # 'DS_GAL9_GLM','DS_RAL3p1_CTC','DS_AfricaLAM'
 #   preprocess_um_data(timeid=int(sys.argv[1]),model='DS_RAL3p1_CTC')

    # Preprocess apparent heat source. 6G, 60mins, 2 tasks. 2-321 inc. for full
#    preprocess_era5_data(timeid=int(sys.argv[1]))

    # Integrate and merge - use sbatch as quite large (20GB mem)
    merge_and_integrate(model='DS_ERA5')


if __name__ == '__main__':

    main()
