load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;;---------------------------------------- read data ---------------------------------------;;
wbt = new((/42,90,45,180/), float)
f1 = addfile("/public/home/fcai/data-observation/ERA5_daily/humid_heat/WetBulbTemp_365days_daymax_1979_2023_2x2.nc", "r")
wbt(:,0:58,:,:) = (f1->wbt(3:,0:58,{-90:0},:))   ;; DJF, 1982-2023
wbt(:,59:89,:,:) = (f1->wbt(3:,334:364,{-90:0},:))




;;-------------------------------- remove zonal mean (1982-2023) ---------------------------;;
wbt_clm = dim_avg_n_Wrap(wbt(0:41,:,:,:), 0)   ;; 1982-2023
wbt_clm_4D = conform(wbt, wbt_clm, (/1,2,3/))
delete(wbt_clm)
wbt_ano = wbt - wbt_clm_4D
copy_VarCoords(wbt, wbt_ano)
delete(wbt_clm_4D)




;;---------------------------------- standardized (1982-2023) -----------------------------;;
wbt_std = dim_stddev_n_Wrap(wbt(0:41,:,:,:), 0)   ;; 1982-2023
wbt_std_4D = conform(wbt, wbt_std, (/1,2,3/))
wbt_std_4D = where(wbt_std_4D.gt.0.0, wbt_std_4D, 1.0e20)
wbt_std_4D@_FillValue = 1.0e20
delete(wbt_std)
copy_VarCoords(wbt_ano, wbt_std_4D)
delete(wbt)


wbt_stz = wbt_ano/wbt_std_4D
copy_VarCoords(wbt_ano, wbt_stz)
delete([/wbt_ano, wbt_std_4D/])
printVarSummary(wbt_stz)
printMinMax(wbt_stz, 1)




;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/"
system("/bin/rm -f "+output+"ERA5_wbt_stz_DJF_1982_2023.nc")
ncdf=addfile(output+"ERA5_wbt_stz_DJF_1982_2023.nc","c")
ncdf->wbt_stz=wbt_stz


end





