load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin


;;---------------------------------------- read data ---------------------------------------;;
f1 = addfile("1_threshold_1x1/ERA5_tmax_stz_JJA_1982_2023.nc", "r")
tmax_stz = (f1->tmax_stz)
printMinMax(tmax_stz, 1)




tmax_95th := new((/91,360/), float)
do ilat = 0,90
print(" ilat = "+ilat)
do jlon = 0,359
tmax_stz_1D = ndtooned(tmax_stz(:,:,ilat,jlon))
tmax_stz_1D = where(tmax_stz_1D.eq.1.0e20, 0, tmax_stz_1D)
qsort(tmax_stz_1D)

tmax_95th(ilat,jlon) = tmax_stz_1D(42*92-194)  ;; top 5%

end do
end do
printMinMax(tmax_95th, 1)


tmax_95th_4D = conform(tmax_stz, tmax_95th, (/2,3/))
tmax95 = where(tmax_stz.gt.tmax_95th_4D, 1, 0)
print(" Sum =  "+ dim_sum(dim_sum(dim_sum(dim_sum(dim_sum(tmax95))))))






;;----------------------------  to nc file  ----------------------------;;
output="1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_tmax95_JJA_1982_2023.nc")
ncdf=addfile(output+"ERA5_tmax95_JJA_1982_2023.nc","c")
ncdf->tmax95=tmax95







end
















