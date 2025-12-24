load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin


;;---------------------------------------- read data ---------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part13_qs95_qs95_NH_lag5/threshold/ERA5_qs_stz_DJF_1982_2023_lead5.nc", "r")
qs_stz = (f1->qs_stz)
printMinMax(qs_stz, 1)
dim_0 = dimsizes(qs_stz)




qs_95th := new((/45,180/), float)
do ilat = 0,44
print(" ilat = "+ilat)
do jlon = 0,179
qs_stz_1D = ndtooned(qs_stz(:,:,ilat,jlon))
qs_stz_1D = where(qs_stz_1D.eq.1.0e20, 0, qs_stz_1D)
qsort(qs_stz_1D)

qs_95th(ilat,jlon) = qs_stz_1D(42*90-190)  ;; top 5%

end do
end do
printMinMax(qs_95th, 1)


qs_95th_4D = conform(qs_stz, qs_95th, (/2,3/))
qs95 = where(qs_stz.gt.qs_95th_4D, 1, 0)
print(" % =  "+ int2flt(dim_sum(dim_sum(dim_sum(dim_sum(dim_sum(qs95)))))) / int2flt((dim_0(0)*dim_0(1)*dim_0(2)*dim_0(3))))






;;----------------------------  to nc file  ----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part13_qs95_qs95_NH_lag5/threshold/"
system("/bin/rm -f "+output+"ERA5_qs95_DJF_1982_2023_lead5.nc")
ncdf=addfile(output+"ERA5_qs95_DJF_1982_2023_lead5.nc","c")
ncdf->qs95=qs95
ncdf->qs_95th=qs_95th







end
















