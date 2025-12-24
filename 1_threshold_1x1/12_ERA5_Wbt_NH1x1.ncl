load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin

;;---------------------------------------- read data ---------------------------------------;;
f1 = addfile("1_threshold_1x1/ERA5_wbt_stz_JJA_1982_2023.nc", "r")
wbt_stz = (f1->wbt_stz)
printMinMax(wbt_stz, 1)




wbt_95th := new((/91,360/), float)
do ilat = 0,90
print(" ilat = "+ilat)
do jlon = 0,359
wbt_stz_1D = ndtooned(wbt_stz(:,:,ilat,jlon))
wbt_stz_1D = where(wbt_stz_1D.eq.1.0e20, 0, wbt_stz_1D)
qsort(wbt_stz_1D)

wbt_95th(ilat,jlon) = wbt_stz_1D(42*92-194)  ;; top 5%

end do
end do
printMinMax(wbt_95th, 1)


wbt_95th_4D = conform(wbt_stz, wbt_95th, (/2,3/))
wbt95 = where(wbt_stz.gt.wbt_95th_4D, 1, 0)
print(" Sum =  "+ dim_sum(dim_sum(dim_sum(dim_sum(dim_sum(wbt95))))))






;;----------------------------  to nc file  ----------------------------;;
output="1_threshold_1x1/"
system("/bin/rm -f "+output+"ERA5_wbt95_JJA_1982_2023.nc")
ncdf=addfile(output+"ERA5_wbt95_JJA_1982_2023.nc","c")
ncdf->wbt95=wbt95







end
















