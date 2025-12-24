load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;------------------------ read lat, lon -------------------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo0 = f00->topo({-90:0},:)
lat = topo0&lat
lon = topo0&lon
dim_0 = dimsizes(topo0)



;;-----------------------------------------------------------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/ERA5_wbt95_DJF_1982_2023.nc","r")
wbt_yearly = int2flt(dim_sum_n_Wrap(f1->wbt95(:,:,:,:), 1))
printMinMax(wbt_yearly, 1)
wbt_yearly@_FillValue = 1.0e20

f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part31_wbt95_sst95_NH/threshold/OISST_sst95_DJF_1982_2023.nc","r")
sst_yearly = int2flt(dim_sum_n_Wrap(f1->sst95(:,:,:,:), 1))
sst_yearly@_FillValue = 1.0e20

wbt_yearly!1 = "lat"
wbt_yearly!2 = "lon"
wbt_yearly&lat = lat
wbt_yearly&lon = lon
copy_VarCoords(wbt_yearly, sst_yearly)



;;--------------------------------------  read distance and networks -----------------------------------------;;
f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_SH.nc","r")
distance = f0->distance(:,:,:,:)


f5 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part32_wbt95_sst95_SH/wbt95_sst95_SH2x2_1982_2023_sign99.nc","r")
networks0 = f5->networks0


;;-------------- most-connected SST Time series (<2,000km) -------------;;
Teleconnection_SST = new((/42,dim_0(0),dim_0(1)/), float)
do ilat = 0,dim_0(0)-1
print(" ilat = "+ilat)
do jlon = 0,dim_0(1)-1
; print(" jlon = "+jlon)

if(.not.ismissing(dim_avg( wbt_yearly(:,ilat,jlon) )))then

networks0_new = where(distance(ilat,jlon,:,:).le.3000.0, networks0(ilat,jlon,:,:), 0.0)
networks0_new_1D = ndtooned(networks0_new)
indices  = ind_resolve(maxind(networks0_new_1D), dim_0)
Teleconnection_SST(:,ilat,jlon) = sst_yearly(:,indices(0,0),indices(0,1))

delete(([/networks0_new,networks0_new_1D,indices/]))
end if
end do
end do
delete(networks0)


copy_VarCoords(wbt_yearly, Teleconnection_SST)






;;------------  to degree file  ------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig3_grid_based_SST/"
name1="Teleconnection_SST_extremes_mostconnected_SH"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->Teleconnection_SST=Teleconnection_SST





end






