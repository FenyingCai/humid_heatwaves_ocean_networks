load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin


;;---------------------------------------- read data ---------------------------------------;;
f1 = addfile("1_threshold_1x1/ERA5_qsurface_stz_DJF_1982_2023.nc", "r")
q_stz = (f1->q_stz)
printMinMax(q_stz, 1)




q_95th := new((/91,360/), float)
do ilat = 0,90
print(" ilat = "+ilat)
do jlon = 0,359
q_stz_1D = ndtooned(q_stz(:,:,ilat,jlon))
q_stz_1D = where(q_stz_1D.eq.1.0e20, 0, q_stz_1D)
qsort(q_stz_1D)

q_95th(ilat,jlon) = q_stz_1D(42*90-190)  ;; top 5%

end do
end do
printMinMax(q_95th, 1)


q_95th_4D = conform(q_stz, q_95th, (/2,3/))
q95 = where(q_stz.gt.q_95th_4D, 1, 0)
print(" Sum =  "+ dim_sum(dim_sum(dim_sum(dim_sum(dim_sum(q95))))))






;;----------------------------  to nc file  ----------------------------;;
output="1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_qsurface95_DJF_1982_2023.nc")
ncdf=addfile(output+"ERA5_qsurface95_DJF_1982_2023.nc","c")
ncdf->q95=q95







end
















