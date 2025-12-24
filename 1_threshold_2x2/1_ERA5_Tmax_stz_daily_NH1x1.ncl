load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin



;;---------------------------------------- read data ---------------------------------------;;
tmax = new((/45,92,91,360/), float)
do i = 0,44
print(" year = " + (i+1979))

f1 = addfile("data/t2m_max/t2m."+flt2string(i+1979)+"-06.daily.nc", "r")
tmax(i,0:29,:,:) = short2flt(f1->t2m(:,{0:90},:))
f2 = addfile("data/t2m_max/t2m."+flt2string(i+1979)+"-07.daily.nc", "r")
tmax(i,30:60,:,:) = short2flt(f2->t2m(:,{0:90},:))
f3 = addfile("data/t2m_max/t2m."+flt2string(i+1979)+"-08.daily.nc", "r")
tmax(i,61:91,:,:) = short2flt(f3->t2m(:,{0:90},:))
end do
printMinMax(tmax, 1)



;;-------------------------------- remove zonal mean (1979-2023) ---------------------------;;
tmax_clm = dim_avg_n_Wrap(tmax(0:44,:,:,:), 0)   ;; 1979-2023
tmax_clm_4D = conform(tmax, tmax_clm, (/1,2,3/))
delete(tmax_clm)
tmax_ano = tmax - tmax_clm_4D
copy_VarCoords(tmax, tmax_ano)
delete(tmax_clm_4D)




;;---------------------------------- standardized (1979-2023) -----------------------------;;
tmax_std = dim_stddev_n_Wrap(tmax(0:44,:,:,:), 0)   ;; 1979-2023
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
system("/bin/rm -f "+output+"ERA5_tmax_stz_JJA_1979_2023.nc")
ncdf=addfile(output+"ERA5_tmax_stz_JJA_1979_2023.nc","c")
ncdf->tmax_stz=tmax_stz


end





