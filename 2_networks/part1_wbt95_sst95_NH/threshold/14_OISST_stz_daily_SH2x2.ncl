load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;;---------------------------------------- read data ---------------------------------------;;
sst = new((/42,90,45,180/), float)
do iyear = 0,41
print("iyear = "+(iyear+1982))
f1 = addfile("/public/home/fcai/data0/OISST_daily2x2/OISST_"+flt2string(iyear+1982)+".nc", "r")
sst(iyear,0:58,:,:) = (f1->sst(0:58,{-90:0},:))   ;; DJF, 1982-2023
sst(iyear,59:89,:,:) = (f1->sst(334:364,{-90:0},:))   ;; DJF, 1982-2023
end do




;;-------------------------------- remove zonal mean (1982-2023) ---------------------------;;
sst_clm = dim_avg_n_Wrap(sst(0:41,:,:,:), 0)   ;; 1982-2023
sst_clm_4D = conform(sst, sst_clm, (/1,2,3/))
delete(sst_clm)
sst_ano = sst - sst_clm_4D
copy_VarCoords(sst, sst_ano)
delete(sst_clm_4D)




;;---------------------------------- standardized (1982-2023) -----------------------------;;
sst_std = dim_stddev_n_Wrap(sst(0:41,:,:,:), 0)   ;; 1982-2023
sst_std_4D = conform(sst, sst_std, (/1,2,3/))
sst_std_4D = where(sst_std_4D.gt.0.0, sst_std_4D, 1.0e20)
sst_std_4D@_FillValue = 1.0e20
delete(sst_std)
copy_VarCoords(sst_ano, sst_std_4D)
delete(sst)


sst_stz = sst_ano/sst_std_4D
copy_VarCoords(sst_ano, sst_stz)
sst_stz!0 = "year"
sst_stz!1 = "day"
delete([/sst_ano, sst_std_4D/])
printVarSummary(sst_stz)
printMinMax(sst_stz, 1)




;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/"
system("/bin/rm -f "+output+"OISST_sst_stz_DJF_1982_2023.nc")
ncdf=addfile(output+"OISST_sst_stz_DJF_1982_2023.nc","c")
ncdf->sst_stz=sst_stz


end





