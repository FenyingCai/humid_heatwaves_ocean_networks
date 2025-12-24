load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;;---------------------------------------- read data ---------------------------------------;;
tmax = new((/42,90,91,360/), float)
do i = 0,41
print(" year = " + (i+1982))

f1 = addfile("data/t2m_max/t2m."+flt2string(i+1982)+"-01.daily.nc", "r")
tmax(i,0:30,:,:) = short2flt(f1->t2m(:,{-90:0},:))
f2 = addfile("data/t2m_max/t2m."+flt2string(i+1982)+"-02.daily.nc", "r")
tmax(i,31:58,:,:) = short2flt(f2->t2m(0:27,{-90:0},:))
f3 = addfile("data/t2m_max/t2m."+flt2string(i+1982)+"-12.daily.nc", "r")
tmax(i,59:89,:,:) = short2flt(f3->t2m(:,{-90:0},:))
end do
printMinMax(tmax, 1)



;;-------------------------------- remove zonal mean (1982-2023) ---------------------------;;
tmax_clm = dim_avg_n_Wrap(tmax(0:41,:,:,:), 0)   ;; 1982-2023
tmax_clm_4D = conform(tmax, tmax_clm, (/1,2,3/))
delete(tmax_clm)
tmax_ano = tmax - tmax_clm_4D
copy_VarCoords(tmax, tmax_ano)
delete(tmax_clm_4D)




;;---------------------------------- standardized (1982-2023) -----------------------------;;
tmax_std = dim_stddev_n_Wrap(tmax(0:41,:,:,:), 0)   ;; 1982-2023
tmax_std_4D = conform(tmax, tmax_std, (/1,2,3/))
tmax_std_4D = where(tmax_std_4D.gt.0.0, tmax_std_4D, 1.0e20)
tmax_std_4D@_FillValue = 1.0e20
delete(tmax_std)
copy_VarCoords(tmax_ano, tmax_std_4D)
delete(tmax)




tmax_stz = tmax_ano/tmax_std_4D
copy_VarCoords(tmax_ano, tmax_stz)
delete([/tmax_ano, tmax_std_4D/])
printVarSummary(tmax_stz)
printMinMax(tmax_stz, 1)




;;--------------------------- to nc file -----------------------------;;
output="1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_tmax_stz_DJF_1982_2023.nc")
ncdf=addfile(output+"ERA5_tmax_stz_DJF_1982_2023.nc","c")
ncdf->tmax_stz=tmax_stz


end





