load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;------------------------------------------------    read lat, lon    --------------------------------------------------;;
; f00 = addfile("/public/home/fcai/data0/topography/topo_5x5.nc","r")
; topo0 = f00->topo({0:80},{0:360})
; lat = topo0&lat
; lon = topo0&lon
; dim_0 = dimsizes(topo0)







;;---------------------------------------- read data ---------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/ERA5_wbt_stz_JJA_1982_2023.nc", "r")
wbt_stz = (f1->wbt_stz)
printMinMax(wbt_stz, 1)
dim_0 = dimsizes(wbt_stz)




wbt_95th := new((/dim_0(2),dim_0(3)/), float)
do ilat = 0,dim_0(2)-1
print(" ilat = "+ilat)
do jlon = 0,dim_0(3)-1
wbt_stz_1D = ndtooned(wbt_stz(:,:,ilat,jlon))
wbt_stz_1D = where(wbt_stz_1D.eq.1.0e20, 0, wbt_stz_1D)
qsort(wbt_stz_1D)

wbt_95th(ilat,jlon) = wbt_stz_1D(42*92-194)  ;; top 5%

end do
end do
printMinMax(wbt_95th, 1)


wbt_95th_4D = conform(wbt_stz, wbt_95th, (/2,3/))
wbt95 = where(wbt_stz.gt.wbt_95th_4D, 1, 0)
print(" Extreme =  "+ dim_sum(dim_sum(dim_sum(dim_sum(wbt95)))))
print(" Sum =  "+ dim_0(0)*dim_0(1)*dim_0(2)*dim_0(3))
print(" % =  "+ int2flt(dim_sum(dim_sum(dim_sum(dim_sum(wbt95)))))/ int2flt(dim_0(0)*dim_0(1)*dim_0(2)*dim_0(3)) )






;;----------------------------  to nc file  ----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/"
system("/bin/rm -f "+output+"ERA5_wbt95_JJA_1982_2023.nc")
ncdf=addfile(output+"ERA5_wbt95_JJA_1982_2023.nc","c")
ncdf->wbt95=wbt95







end
















