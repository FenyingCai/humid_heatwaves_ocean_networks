load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;-----------------------------------------------------------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({-90:90},:)
dim_00 := dimsizes(topo)
degree = new((/42,dim_00(0),dim_00(1)/), float)



do iyear = 0,41
print("iyear = "+iyear)
;;----------------------------------------------------  read networks (NH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part3_wbt95_qs95_NH/wbt95_qs95_NH2x2_1982_2023_sign99.nc","r")
networks_42yrs := (f1->networks0)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part3_wbt95_qs95_NH/wbt95_qs95_NH2x2_"+int2flt(iyear+1982)+".nc","r")
networks1 := (f1->networks0)
networks1 = where(networks_42yrs.gt.0.0, networks1, 0.0)
delete(networks_42yrs)
dim_0 := dimsizes(networks1)


;;------------- remove ocean ---------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo0 = f00->topo({0:90},:)
topo_networks1 = conform(networks1, topo0, (/0,1/))
topo_networks2 = conform(networks1, topo0, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
; printMinMax(networks1, 0)


f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_NH.nc","r")
distance = f0->distance(:,:,:,:)
networks1 = where(distance.le.2000.0, networks1, 1.0e20)
delete(distance)

;;------------------------  lat weighted  ------------------------;;
lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_matrix(i,:,j,:) = cos((i*2+1+j*2+1) /2.0 /180.0 * 3.14159)
  end do
end do

networks1__ = where(networks1.gt.0.0, 1.0, 0.0)
ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 1), 0) 
ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 3), 2) 
ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
ref1@_FillValue = 1.0e20
ref2@_FillValue = 1.0e20
degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2
delete(networks1)
degree_NH = where(topo0.ge.0.0, degree_land, 0.0)







;;----------------------------------------------------  read networks (SH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part4_wbt95_qs95_SH/wbt95_qs95_SH2x2_1982_2023_sign99.nc","r")
networks_42yrs := (f1->networks0)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part4_wbt95_qs95_SH/wbt95_qs95_SH2x2_"+int2flt(iyear+1982)+".nc","r")
networks1 := (f1->networks0)
networks1(0:14,:,:,:) = 0.0
networks1(:,:,0:14,:) = 0.0
networks1 = where(networks_42yrs.gt.0.0, networks1, 0.0)
delete(networks_42yrs)
dim_0 := dimsizes(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo0 = f00->topo({-90:0},:)
topo_networks1 = conform(networks1, topo0, (/0,1/))
topo_networks2 = conform(networks1, topo0, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
; printMinMax(networks1, 0)


f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_SH.nc","r")
distance = f0->distance(:,:,:,:)
networks1 = where(distance.le.2000.0, networks1, 1.0e20)
delete(distance)

;;------------------------  lat weighted  ------------------------;;
lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_matrix(i,:,j,:) = cos((i*2-89+j*2-89) /2.0 /180.0 * 3.14159)
  end do
end do

networks1__ = where(networks1.gt.0.0, 1.0, 0.0)
ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 1), 0) 
ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 3), 2) 
ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
ref1@_FillValue = 1.0e20
ref2@_FillValue = 1.0e20
degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2

degree_SH = where(topo0.ge.0.0, degree_land, 0.0)

;;---------------------  Degree (NH + SH)  ----------------------;;
degree(iyear,0:44,:) = degree_SH * 1.0
degree(iyear,45:89,:) = degree_NH * 1.0
copy_VarCoords(topo, degree(iyear,:,:))
printMinMax(degree(iyear,:,:),1)
end do







printVarSummary(degree)
printMinMax(degree, 1)


;;---------------------------------- to nc file ------------------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_local_concurrence/"
name1="wbt95_qs95_LCle2000_42yrs"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->degree=degree












end







