import numpy as np
import xarray as xr


ds1 = xr.open_dataset("/public/home/fcai/data0/topography/topo_2x2.nc")
topo0 = ds1.topo.loc[-90:0,:]
lat = topo0.lat; lon = topo0.lon



##--##--##--##--                     --##--##--##--##
##--##--##--##--  1982-2023, 42 yrs  --##--##--##--##
##--##--##--##--                     --##--##--##--##

for iyear in range(0,42):
  print("                   ")
  print('iyear = ',iyear+1982)
  print("                   ")
  


  ##--##--##--##--  read 90th threshold   --##--##--##--##
  ds1 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part3_wbt95_qs95_NH/threshold/ERA5_qs95_DJF_1982_2023.nc")
  qs95 = ds1.qs95[iyear,:,:,:]
  

  ds2 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part13_qs95_qs95_NH_lag5/threshold/ERA5_qs95_DJF_1982_2023_lead5.nc")
  qs95_lead5 = ds2.qs95[iyear,:,:,:]




  ##--##--##--##--  network  --##--##--##--##
  hw_networks0 = np.zeros((np.shape(qs95)[1],np.shape(qs95)[2], np.shape(qs95)[1],np.shape(qs95)[2]), np.float32)
  # print(hw_networks0.size * hw_networks0.itemsize / 1073741824.0)

  for alat in range(np.shape(qs95)[1]):
    print('alat = ',alat)
    for alon in range(np.shape(qs95)[2]):

      t_True_3D = np.tile(qs95[:,alat,alon], (np.shape(qs95)[1],np.shape(qs95)[2],1)).transpose(2,0,1)   ## day|92* lat|45* lon|180
      
      concurrent_days = np.sum(t_True_3D[:,:,:]+qs95_lead5[:,:,:]==2,axis=0)
      hw_networks0[:,:,alat,alon] = concurrent_days

      del  t_True_3D, concurrent_days
  del qs95,qs95_lead5
  print("network ", np.max(hw_networks0), np.min(hw_networks0))




  ##--##--##--##--  to nc file  --##--##--##--##
  hw_networks_array0= xr.DataArray(data=hw_networks0.data, dims=['lat1', 'lon1', 'lat2', 'lon2'],
                                  coords={'lat1':lat.data,'lon1':lon.data,
                                          'lat2':lat.data,'lon2':lon.data})
  ds0 = xr.Dataset(data_vars=dict(networks0=hw_networks_array0))
  ds0.to_netcdf("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part14_qs95_qs95_SH_lag5/qs95_qs95_lag5_SH2x2_"+str(iyear+1982)+".nc")
  ds0.close()
  del hw_networks0, hw_networks_array0




