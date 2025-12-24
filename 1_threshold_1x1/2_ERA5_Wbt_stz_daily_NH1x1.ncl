load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;;---------------------------------------- read data ---------------------------------------;;
wbt = new((/42,92,91,360/), float)
f1 = addfile("1_threshold_1x1/WetBulbTemp_365days_daymax_1979_2023.nc", "r")
wbt(:,:,:,:) = (f1->wbt(3:44,151:151+91,90:180,:))   ;; JJA, 0-90N
lat=fspan(0.0, 90.0, 91)
lat@units="degrees_north"
wbt&latitude = lat



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
output="1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_wbt_stz_JJA_1982_2023.nc")
ncdf=addfile(output+"ERA5_wbt_stz_JJA_1982_2023.nc","c")
ncdf->wbt_stz=wbt_stz


end





