load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin






;;------------------------------------------------- NH ------------------------------------------------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({0:90},:)
dim_0 = dimsizes(topo)

latt=topo&lat
lonn=topo&lon
;;--------------------------  lat weighted  ---------------------------;;
lat_weight = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_weight(i,:,j,:) = cos((i*2+1+j*2+1) /2.0 /180.0 * 3.14159)
  end do
end do
printMinMax(lat_weight, 1)


;;---------------------------------  read three oceans  -----------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/Distance_landocean_NH.nc","r")
land1 = f1->land1
land2 = f1->land2
ocean1 = f1->ocean1
ocean2 = f1->ocean2

lat_weight2D = new((/dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  lat_weight2D(i,:) = cos((i*2+1) /180.0 * 3.14159)
end do
area_land1 = dim_sum_n_Wrap(dim_sum_n_Wrap(land1*lat_weight2D, 1), 0)
area_land2 = dim_sum_n_Wrap(dim_sum_n_Wrap(land2*lat_weight2D, 1), 0)
area_ocean1 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean1*lat_weight2D, 1), 0)
area_ocean2 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean2*lat_weight2D, 1), 0)
print("area = "+area_ocean2+"  "+area_ocean1+"  "+area_land1+"  "+area_land2)


land1_from = conform(lat_weight, land1, (/0,1/))
land1_to = conform(lat_weight, land1, (/2,3/))
land2_from = conform(lat_weight, land2, (/0,1/))
land2_to = conform(lat_weight, land2, (/2,3/))
ocean1_from = conform(lat_weight, ocean1, (/0,1/))
ocean1_to = conform(lat_weight, ocean1, (/2,3/))
ocean2_from = conform(lat_weight, ocean2, (/0,1/))
ocean2_to = conform(lat_weight, ocean2, (/2,3/))


;;-----------------------------------  read networks (0-36 months)  -----------------------------------;;
flow = new((/2,4,4/), float)
networks = new((/2,dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part11_qs95_qs95_NH/qs95_qs95_NH2x2_1982_2023_sign99.nc","r")
networks(0,:,:,:,:) = (f1->networks0)
f6 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part13_qs95_qs95_NH_lag5/qs95_qs95_lag5_NH2x2_1982_2023_sign99.nc","r")
networks(1,:,:,:,:) = (f6->networks0)


do ilag = 0,1
print(" ilag  =  "+ilag)

networks0 := networks(ilag,:,:,:,:)
networks0 = where(ismissing(networks0), 0.0, networks0)
printMinMax(networks0, 1)


flow(ilag,0,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,0,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,0,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*land1_to, 3), 2), 1), 0)
flow(ilag,0,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*land2_to, 3), 2), 1), 0)

flow(ilag,1,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,1,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,1,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*land1_to, 3), 2), 1), 0)
flow(ilag,1,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*land2_to, 3), 2), 1), 0)

flow(ilag,2,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,2,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,2,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*land1_to, 3), 2), 1), 0)
flow(ilag,2,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*land2_to, 3), 2), 1), 0)

flow(ilag,3,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,3,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,3,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*land1_to, 3), 2), 1), 0)
flow(ilag,3,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*land2_to, 3), 2), 1), 0)

printMinMax(flow(ilag,:,:), 1)

end do

flow_weight = flow
flow_weight(:,0,:) = flow_weight(:,0,:) / area_ocean2
flow_weight(:,1,:) = flow_weight(:,1,:) / area_ocean1
flow_weight(:,2,:) = flow_weight(:,2,:) / area_land1
flow_weight(:,3,:) = flow_weight(:,3,:) / area_land2


;;------------------------------ to nc file --------------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/"
name1="Flow_landocean_qs_qs_lag0_5_NH"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->flow=flow
ncdf1->flow_weight=flow_weight












;;------------------------------------------------- SH ------------------------------------------------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:0},:)
dim_0 = dimsizes(topo)

latt=topo&lat
lonn=topo&lon
;;--------------------------  lat weighted  ---------------------------;;
lat_weight = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_weight(i,:,j,:) = cos((i*2-89+j*2-89) /2.0 /180.0 * 3.14159)
  end do
end do
printMinMax(lat_weight, 1)


;;---------------------------------  read three oceans  -----------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/Distance_landocean_SH.nc","r")
land1 = f1->land1
land2 = f1->land2
ocean1 = f1->ocean1
ocean2 = f1->ocean2

lat_weight2D = new((/dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  lat_weight2D(i,:) = cos((i*2-89) /180.0 * 3.14159)
end do
area_land1 = dim_sum_n_Wrap(dim_sum_n_Wrap(land1*lat_weight2D, 1), 0)
area_land2 = dim_sum_n_Wrap(dim_sum_n_Wrap(land2*lat_weight2D, 1), 0)
area_ocean1 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean1*lat_weight2D, 1), 0)
area_ocean2 = dim_sum_n_Wrap(dim_sum_n_Wrap(ocean2*lat_weight2D, 1), 0)
print("area = "+area_ocean2+"  "+area_ocean1+"  "+area_land1+"  "+area_land2)


land1_from = conform(lat_weight, land1, (/0,1/))
land1_to = conform(lat_weight, land1, (/2,3/))
land2_from = conform(lat_weight, land2, (/0,1/))
land2_to = conform(lat_weight, land2, (/2,3/))
ocean1_from = conform(lat_weight, ocean1, (/0,1/))
ocean1_to = conform(lat_weight, ocean1, (/2,3/))
ocean2_from = conform(lat_weight, ocean2, (/0,1/))
ocean2_to = conform(lat_weight, ocean2, (/2,3/))


;;-----------------------------------  read networks (0-36 months)  -----------------------------------;;
flow = new((/2,4,4/), float)
networks = new((/2,dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part12_qs95_qs95_SH/qs95_qs95_SH2x2_1982_2023_sign99.nc","r")
networks(0,:,:,:,:) = (f1->networks0)
f6 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part14_qs95_qs95_SH_lag5/qs95_qs95_lag5_SH2x2_1982_2023_sign99.nc","r")
networks(1,:,:,:,:) = (f6->networks0)


do ilag = 0,1
print(" ilag  =  "+ilag)

networks0 := networks(ilag,:,:,:,:)
networks0 = where(ismissing(networks0), 0.0, networks0)
printMinMax(networks0, 1)



flow(ilag,0,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,0,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,0,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*land1_to, 3), 2), 1), 0)
flow(ilag,0,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean2_from*land2_to, 3), 2), 1), 0)

flow(ilag,1,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,1,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,1,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*land1_to, 3), 2), 1), 0)
flow(ilag,1,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * ocean1_from*land2_to, 3), 2), 1), 0)

flow(ilag,2,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,2,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,2,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*land1_to, 3), 2), 1), 0)
flow(ilag,2,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land1_from*land2_to, 3), 2), 1), 0)

flow(ilag,3,0) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*ocean2_to, 3), 2), 1), 0)
flow(ilag,3,1) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*ocean1_to, 3), 2), 1), 0)
flow(ilag,3,2) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*land1_to, 3), 2), 1), 0)
flow(ilag,3,3) = dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(dim_sum_n_Wrap(networks0*lat_weight * land2_from*land2_to, 3), 2), 1), 0)

printMinMax(flow(ilag,:,:), 1)

end do

flow_weight = flow
flow_weight(:,0,:) = flow_weight(:,0,:) / area_ocean2
flow_weight(:,1,:) = flow_weight(:,1,:) / area_ocean1
flow_weight(:,2,:) = flow_weight(:,2,:) / area_land1
flow_weight(:,3,:) = flow_weight(:,3,:) / area_land2


;;------------------------------ to nc file --------------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/"
name1="Flow_landocean_qs_qs_lag0_5_SH"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->flow=flow
ncdf1->flow_weight=flow_weight








end


