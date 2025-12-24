load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;------------------------------------------------ read topo ------------------------------------------------;;
f1 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc", "r")
topo0 = (f1->topo({-90:0},:))
dim_0 = dimsizes(topo0)



;;---------------------------------------- read (qs) ---------------------------------------;;
qs = new((/42,90,45,180/), float)
do i = 0,41
print(" year = " + (i+1982))
f2 = addfile("/public/home/fcai/data0/ERA5_q2m_daily2x2/q2m."+flt2string(i+1981)+".nc", "r")
qs(i,0:4,:,:) = (f2->q2m(360:364,{-90:0},:))

f1 = addfile("/public/home/fcai/data0/ERA5_q2m_daily2x2/q2m."+flt2string(i+1982)+".nc", "r")
qs(i,0+5:58,:,:) = (f1->q2m(0:58-5,{-90:0},:))   ;; DJF, lead = 5
qs(i,59:89,:,:) = (f1->q2m(334-5:364-5,{-90:0},:))
end do





;;-------------------------------- remove temporal mean (1982-2023) ---------------------------;;
qs_clm = dim_avg_n_Wrap(qs(0:41,:,:,:), 0)   ;; 1982-2023
qs_clm_4D = conform(qs, qs_clm, (/1,2,3/))
delete(qs_clm)
qs_ano = qs - qs_clm_4D
copy_VarCoords(qs, qs_ano)
delete(qs_clm_4D)




;;---------------------------------- standardized (1982-2023) -----------------------------;;
qs_std = dim_stddev_n_Wrap(qs(0:41,:,:,:), 0)   ;; 1982-2023
qs_std_4D = conform(qs, qs_std, (/1,2,3/))
qs_std_4D = where(qs_std_4D.gt.0.0, qs_std_4D, 1.0e20)
qs_std_4D@_FillValue = 1.0e20
delete(qs_std)
copy_VarCoords(qs_ano, qs_std_4D)



qs_stz = qs_ano/qs_std_4D
copy_VarCoords(qs_ano, qs_stz)
delete([/qs_ano, qs_std_4D/])
printVarSummary(qs_stz)
printMinMax(qs_stz, 1)




;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part13_qs95_qs95_NH_lag5/threshold/"
system("/bin/rm -f "+output+"ERA5_qs_stz_DJF_1982_2023_lead5.nc")
ncdf=addfile(output+"ERA5_qs_stz_DJF_1982_2023_lead5.nc","c")
ncdf->qs_stz=qs_stz
ncdf->qs=qs


end





