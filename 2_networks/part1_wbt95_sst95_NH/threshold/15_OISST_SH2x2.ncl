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






; wbt95 = new((/45,92,91,360/), integer)

;;---------------------------------------- read data ---------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/OISST_sst_stz_DJF_1982_2023.nc", "r")
sst_stz = (f1->sst_stz)
printMinMax(sst_stz, 1)
dim_0 = dimsizes(sst_stz)




sst_95th := new((/dim_0(2),dim_0(3)/), float)
do ilat = 0,dim_0(2)-1
print(" ilat = "+ilat)
do jlon = 0,dim_0(3)-1
sst_stz_1D = ndtooned(sst_stz(:,:,ilat,jlon))
sst_stz_1D = where(sst_stz_1D.eq.1.0e20, 0, sst_stz_1D)
qsort(sst_stz_1D)

sst_95th(ilat,jlon) = sst_stz_1D(42*90-190)  ;; top 5%

end do
end do
printMinMax(sst_95th, 1)


sst_95th_4D = conform(sst_stz, sst_95th, (/2,3/))
sst95 = where(sst_stz.gt.sst_95th_4D, 1, 0)
print(" Extreme =  "+ dim_sum(dim_sum(dim_sum(dim_sum(sst95)))))
print(" Sum =  "+ dim_0(0)*dim_0(1)*dim_0(2)*dim_0(3))
print(" % =  "+ int2flt(dim_sum(dim_sum(dim_sum(dim_sum(sst95)))))/ int2flt(dim_0(0)*dim_0(1)*dim_0(2)*dim_0(3)) )
do ilat = 0,dim_0(2)-1,5
do jlon = 0,dim_0(3)-1,5
print(" % =  "+ int2flt(dim_sum(dim_sum(sst95(:,:,ilat,jlon))))/ int2flt(dim_0(0)*dim_0(1)) )
end do
end do





;;----------------------------  to nc file  ----------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/"
system("/bin/rm -f "+output+"OISST_sst95_DJF_1982_2023.nc")
ncdf=addfile(output+"OISST_sst95_DJF_1982_2023.nc","c")
ncdf->sst95=sst95







end
















