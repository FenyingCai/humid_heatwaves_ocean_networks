load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;;---------------------------------------- read data ---------------------------------------;;
wbt = new((/45,90,91,360/), float)
f1 = addfile("/public/home/fcai/data-observation/ERA5_daily/humid_heat/WetBulbTemp_365days_daymax_1979_2023.nc", "r")
wbt(:,0:58,:,:) = (f1->wbt(:,0:58,0:90,:))   ;; DJF, 90S-0
wbt(:,59:89,:,:) = (f1->wbt(:,334:364,0:90,:))
lat=fspan(-90.0, 0.0, 91)
lat@units="degrees_north"
wbt&latitude = lat



;;-------------------------------- remove zonal mean (1979-2023) ---------------------------;;
wbt_clm = dim_avg_n_Wrap(wbt(0:44,:,:,:), 0)   ;; 1979-2023
wbt_clm_4D = conform(wbt, wbt_clm, (/1,2,3/))
delete(wbt_clm)
wbt_ano = wbt - wbt_clm_4D
copy_VarCoords(wbt, wbt_ano)
delete(wbt_clm_4D)




;;---------------------------------- standardized (1979-2023) -----------------------------;;
wbt_std = dim_stddev_n_Wrap(wbt(0:44,:,:,:), 0)   ;; 1979-2023
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
output="/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_wbt_stz_DJF_1979_2023.nc")
ncdf=addfile(output+"ERA5_wbt_stz_DJF_1979_2023.nc","c")
ncdf->wbt_stz=wbt_stz


end





