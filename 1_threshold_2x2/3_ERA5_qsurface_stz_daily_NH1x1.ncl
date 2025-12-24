load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;---------------------------------------- read ps ---------------------------------------;;
ps = new((/45,92,91,360/), float)
do i = 0,44
print(" year = " + (i+1979))

f1 = addfile("/public/home/fcai/data-observation/ERA5_daily/ps/ps."+flt2string(i+1979)+"-06.daily.nc", "r")
ps(i,0:29,:,:) = short2flt(f1->sp(:,{0:90},:))
f2 = addfile("/public/home/fcai/data-observation/ERA5_daily/ps/ps."+flt2string(i+1979)+"-07.daily.nc", "r")
ps(i,30:60,:,:) = short2flt(f2->sp(:,{0:90},:))
f3 = addfile("/public/home/fcai/data-observation/ERA5_daily/ps/ps."+flt2string(i+1979)+"-08.daily.nc", "r")
ps(i,61:91,:,:) = short2flt(f3->sp(:,{0:90},:))
end do

ps_clm = dim_avg_n_Wrap(dim_avg_n_Wrap(ps, 1), 0) / 100.0
printMinMax(ps_clm, 0)
; print(ps_clm)

lev = (/1000,925,850,700,600,500/)
lev_3D = conform(ps(0,0:5,:,:), lev, (/0/))
ps_3D = conform(ps(0,0:5,:,:), ps_clm, (/1,2/))
delete(ps)
ps_diff = (ps_3D - lev_3D)
ps_diff = where(ps_diff.ge.0.0, ps_diff, 1.0e20)
ps_diff@_FillValue = 1.0e20

ps_index = dim_minind(ps_diff, 0)
printMinMax(ps_index, 1)





;;--------------------------------------------- read data ---------------------------------------------;;
q0 = new((/45,92,6,91,360/), float)
do i = 0,44
print(" year = " + (i+1979))

f1 = addfile("/public/home/fcai/data-observation/ERA5_daily/q/q."+flt2string(i+1979)+"-06.daily.nc", "r")
q0(i,0:29,:,:,:) = short2flt(f1->q(:,{lev},{0:90},:))
f2 = addfile("/public/home/fcai/data-observation/ERA5_daily/q/q."+flt2string(i+1979)+"-07.daily.nc", "r")
q0(i,30:60,:,:,:) = short2flt(f2->q(:,{lev},{0:90},:))
f3 = addfile("/public/home/fcai/data-observation/ERA5_daily/q/q."+flt2string(i+1979)+"-08.daily.nc", "r")
q0(i,61:91,:,:,:) = short2flt(f3->q(:,{lev},{0:90},:))
end do




q = new((/45,92,91,360/), float)
do ilat = 0,90
print(" ilat = "+ilat)
do ilon = 0,359

q(:,:,ilat,ilon) = q0(:,:,ps_index(ilat,ilon),ilat,ilon)

end do
end do
delete(q0)




;;-------------------------------- remove zonal mean (1979-2023) ---------------------------;;
q_clm = dim_avg_n_Wrap(q(0:44,:,:,:), 0)   ;; 1979-2023
q_clm_4D = conform(q, q_clm, (/1,2,3/))
delete(q_clm)
q_ano = q - q_clm_4D
copy_VarCoords(q, q_ano)
delete(q_clm_4D)




;;---------------------------------- standardized (1979-2023) -----------------------------;;
q_std = dim_stddev_n_Wrap(q(0:44,:,:,:), 0)   ;; 1979-2023
q_std_4D = conform(q, q_std, (/1,2,3/))
q_std_4D = where(q_std_4D.gt.0.0, q_std_4D, 1.0e20)
q_std_4D@_FillValue = 1.0e20
delete(q_std)
copy_VarCoords(q_ano, q_std_4D)
delete(q)


q_stz = q_ano/q_std_4D
copy_VarCoords(q_ano, q_stz)
delete([/q_ano, q_std_4D/])
printVarSummary(q_stz)
printMinMax(q_stz, 1)




;;--------------------------- to nc file -----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_qsurface_stz_JJA_1979_2023.nc")
ncdf=addfile(output+"ERA5_qsurface_stz_JJA_1979_2023.nc","c")
ncdf->q_stz=q_stz


end





