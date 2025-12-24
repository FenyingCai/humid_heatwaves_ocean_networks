load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;;--------------------------------- read data ---------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 = (f1->topo({-90:90},:))

; lat=fspan(1,89.0,45)
; lon=fspan(1,359.0,180)
; lat@units="degrees_north"
; lon@units="degrees_east"



f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/figz_similarity/ERA5_wbt95_clm_1982_2023.nc", "r")
wbt95_clm = (f3->wbt95_clm)
copy_VarCoords(topo0, wbt95_clm)








;;----------------------------------------------------  read networks (NH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part7_wbt95_wbt95_NH/wbt95_wbt95_NH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1 = where(networks1.gt.0.0, networks1, 0.0)
printVarSummary(networks1)
printMinMax(networks1, 1)
dim_0 = dimsizes(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({0:90},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
printMinMax(networks1, 0)


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
degree_NH = where(topo.ge.0.0, degree_land, 0.0)







;;----------------------------------------------------  read networks (SH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part8_wbt95_wbt95_SH/wbt95_wbt95_SH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1(0:14,:,:,:) = 0.0
networks1(:,:,0:14,:) = 0.0
networks1 = where(networks1.gt.0.0, networks1, 0.0)
printVarSummary(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:0},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
printMinMax(networks1, 0)


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

degree_SH = where(topo.ge.0.0, degree_land, 0.0)

;;---------------------  Degree (NH + SH)  ----------------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({-90:90},:)
dim_0 := dimsizes(topo)
degree = new((/dim_0(0),dim_0(1)/), float)
degree(0:44,:) = degree_SH * 1.0
degree(45:89,:) = degree_NH * 1.0
copy_VarCoords(topo, degree)
printMinMax(degree,1)


















;;---------------------------------------------------------------------------;;
degree_1D = ndtooned(degree)
degree_1D@_FillValue = -999
qsort(degree_1D)
print(dimsizes(degree_1D)) ;; 16200
degree_1D_ = where(degree_1D.gt.0.0, 1.0, 0.0)
print(dim_sum(degree_1D_)) ;; 3849
print(degree_1D(16200-1-1283))


wbt95_clm_1D = ndtooned(wbt95_clm)
wbt95_clm_1D@_FillValue = -999
print(dimsizes(wbt95_clm_1D)) ;; 65160
wbt95_clm_1D_ = where(wbt95_clm_1D.gt.0.0, 1.0, 0.0)
print(dim_sum(wbt95_clm_1D_)) ;; 15032
qsort(wbt95_clm_1D)
print(wbt95_clm_1D(65160-1-3010))



wbt95_clm_2x2 = linint2(wbt95_clm&lon,wbt95_clm&lat, wbt95_clm, True, topo&lon,topo&lat, 0)
wbt95_clm_2x2@_FillValue=1e+20

both_high = where(wbt95_clm_2x2.gt.25.82957.and.degree.gt.28.48318, 1.0, 0.0)
copy_VarCoords(topo, both_high)












; ;;----------------------------------------------------  read networks (NH)  ----------------------------------------------------------;;
; f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part61_qs95_qs95_NH/qs95_qs95_NH2x2_1982_2023_sign99.nc","r")
; networks1 = (f1->networks0)
; networks1 = where(networks1.gt.0.0, networks1, 0.0)
; printVarSummary(networks1)
; printMinMax(networks1, 1)
; dim_0 := dimsizes(networks1)

; ;;------------- remove ocean ---------------;;
; f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
; topo := f00->topo({0:90},:)
; topo_networks1 = conform(networks1, topo, (/0,1/))
; topo_networks2 = conform(networks1, topo, (/2,3/))

; networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
; networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
; network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
; network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
; delete([/topo_networks1,topo_networks2/])
; networks1@_FillValue = 1.0e20
; printMinMax(networks1, 0)


; f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_NH.nc","r")
; distance = f0->distance(:,:,:,:)
; networks1 = where(distance.le.2000.0, networks1, 1.0e20)
; delete(distance)

; ;;------------------------  lat weighted  ------------------------;;
; lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
; do i = 0,dim_0(0)-1
;   do j = 0,dim_0(0)-1
;     lat_matrix(i,:,j,:) = cos((i*2+1+j*2+1) /2.0 /180.0 * 3.14159)
;   end do
; end do

; networks1__ = where(networks1.gt.0.0, 1.0, 0.0)
; ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 1), 0) 
; ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 3), 2) 
; ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
; ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
; ref1@_FillValue = 1.0e20
; ref2@_FillValue = 1.0e20
; degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2
; delete(networks1)
; degree_NH = where(topo.ge.0.0, degree_land, 0.0)


; ;;----------------------------------------------------  read networks (SH)  ----------------------------------------------------------;;
; f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part62_qs95_qs95_SH/qs95_qs95_SH2x2_1982_2023_sign99.nc","r")
; networks1 = (f1->networks0)
; networks1(0:14,:,:,:) = 0.0
; networks1(:,:,0:14,:) = 0.0
; networks1 = where(networks1.gt.0.0, networks1, 0.0)
; printVarSummary(networks1)

; ;;------------- remove ocean ---------------;;
; f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
; topo = f00->topo({-90:0},:)
; topo_networks1 = conform(networks1, topo, (/0,1/))
; topo_networks2 = conform(networks1, topo, (/2,3/))

; networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
; networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
; network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
; network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
; delete([/topo_networks1,topo_networks2/])
; networks1@_FillValue = 1.0e20
; printMinMax(networks1, 0)


; f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_SH.nc","r")
; distance = f0->distance(:,:,:,:)
; networks1 = where(distance.le.2000.0, networks1, 1.0e20)
; delete(distance)

; ;;------------------------  lat weighted  ------------------------;;
; lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
; do i = 0,dim_0(0)-1
;   do j = 0,dim_0(0)-1
;     lat_matrix(i,:,j,:) = cos((i*2-89+j*2-89) /2.0 /180.0 * 3.14159)
;   end do
; end do

; networks1__ = where(networks1.gt.0.0, 1.0, 0.0)
; ref1 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 1), 0) 
; ref2 = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1__(:,:,:,:)*network_ref, 3), 2) 
; ref1 = where(ref1.gt.0.0, ref1, 1.0e20)
; ref2 = where(ref2.gt.0.0, ref2, 1.0e20)
; ref1@_FillValue = 1.0e20
; ref2@_FillValue = 1.0e20
; degree_land = dim_sum_n_Wrap(dim_sum_n_Wrap(networks1(:,:,:,:)*lat_matrix, 3), 2) / ref2

; degree_SH = where(topo.ge.0.0, degree_land, 0.0)

; ;;---------------------  Degree (NH + SH)  ----------------------;;
; f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
; topo := f00->topo({-90:90},:)
; dim_0 := dimsizes(topo)
; degree_Wet_Wet = new((/dim_0(0),dim_0(1)/), float)
; degree_Wet_Wet(0:44,:) = degree_SH * 1.0
; degree_Wet_Wet(45:89,:) = degree_NH * 1.0
; copy_VarCoords(topo, degree_Wet_Wet)
; printMinMax(degree_Wet_Wet,1)













;;----------------------------------------------------  read networks (NH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part3_wbt95_qs95_NH/wbt95_qs95_NH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1 = where(networks1.gt.0.0, networks1, 0.0)
printVarSummary(networks1)
printMinMax(networks1, 1)
dim_0 := dimsizes(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({0:90},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
printMinMax(networks1, 0)


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
degree_NH = where(topo.ge.0.0, degree_land, 0.0)


;;----------------------------------------------------  read networks (SH)  ----------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part4_wbt95_qs95_SH/wbt95_qs95_SH2x2_1982_2023_sign99.nc","r")
networks1 = (f1->networks0)
networks1(0:14,:,:,:) = 0.0
networks1(:,:,0:14,:) = 0.0
networks1 = where(networks1.gt.0.0, networks1, 0.0)
printVarSummary(networks1)

;;------------- remove ocean ---------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({-90:0},:)
topo_networks1 = conform(networks1, topo, (/0,1/))
topo_networks2 = conform(networks1, topo, (/2,3/))

networks1 = where(topo_networks1.ge.0.0, networks1, 1.0e20)
networks1 = where(topo_networks2.ge.0.0, networks1, 1.0e20)
network_ref = where(topo_networks1.ge.0.0, 1.0, 1.0e20)
network_ref = where(topo_networks2.ge.0.0, network_ref, 1.0e20)
delete([/topo_networks1,topo_networks2/])
networks1@_FillValue = 1.0e20
printMinMax(networks1, 0)


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

degree_SH = where(topo.ge.0.0, degree_land, 0.0)

;;---------------------  Degree (NH + SH)  ----------------------;;
f00 = addfile("/public/data2/model/fcai/data0/topography/topo_2x2.nc","r")
topo := f00->topo({-90:90},:)
dim_0 := dimsizes(topo)
degree_HWs_Wet = new((/dim_0(0),dim_0(1)/), float)
degree_HWs_Wet(0:44,:) = degree_SH * 1.0
degree_HWs_Wet(45:89,:) = degree_NH * 1.0
copy_VarCoords(topo, degree_HWs_Wet)
printMinMax(degree_HWs_Wet,1)








;;---------------------  Interannual correlations  ----------------------;;
f3 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_local_concurrence/wbt95_wbt95_LCle2000_42yrs.nc","r")
LC_wbt95_wbt95 := f3->degree
f4 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig1_local_concurrence/wbt95_qs95_LCle2000_42yrs.nc","r")
LC_wbt95_qs95 := f4->degree

wbt = new((/42,dim_0(0),dim_0(1)/), float)
f5 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/threshold/ERA5_wbt95_JJA_1982_2023.nc","r")
wbt(:,0:44,:) = dim_sum_n_Wrap(f5->wbt95, 1)
f6 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part1_wbt95_sst95_NH/threshold/ERA5_wbt95_DJF_1982_2023.nc","r")
wbt(:,45:89,:) = dim_sum_n_Wrap(f6->wbt95, 1)
wbt = where(wbt.gt.0.0, wbt, 1.0e20)
wbt@_FillValue = 1.0e20

Cor = escorc_n(LC_wbt95_wbt95/wbt, LC_wbt95_qs95/wbt, 0, 0)
Cor = where(Cor.gt.0.0, Cor, 0.0)
copy_VarCoords(topo, Cor)
printMinMax(Cor, 1)
print(dim_sum(dim_sum(Cor)))

Cor0 = escorc_n(LC_wbt95_wbt95, LC_wbt95_qs95, 0, 0)
Cor0 = where(Cor0.gt.0.0, Cor0, 0.0)
copy_VarCoords(topo, Cor0)
printMinMax(Cor0, 1)
print(dim_sum(dim_sum(Cor0)))




;;---------------------------------------------------  plots  -----------------------------------------------------
res = True
res@tmBorderThicknessF = 0.8
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnAddCyclic =  True
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.04
res@gsnRightStringFontHeightF  = 0.04
res@gsnCenterStringFontHeightF  = 0.04
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
cnres = res
cnres@gsnAddCyclic =  False
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBOn        = True
res@tmYLOn        = True
; res@tmXBLabelFontHeightF = 0.038
; res@tmYLLabelFontHeightF = 0.038
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0





res@mpMinLatF = -60
res@mpMaxLatF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180-25
res@tmXBMode = "Explicit"
res@tmXBValues = (/-135,-90,-45,0,45,90,135,180,225,270,315/) + 25
res@tmXBLabels = (/"135W","90W","45W","0","45E","90E","135E","180","135W","90W","45W"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-360, 360, 49) + 25
res@tmYLMode = "Explicit"
res@tmYLValues = (/-90,-60,-30,0,30,60,90/)
res@tmYLLabels = (/"90S","60S","30S","0","30N","60N","90N"/)
res@tmYLMinorValues = fspan(-90, 90, 19)
res@tmXBLabelsOn = True
res@tmXBLabelFontHeightF = 0.017
res@tmYLLabelsOn = True
res@tmYLLabelFontHeightF = 0.017

res@tmYLMajorThicknessF = 1.1
res@tmXBMajorThicknessF = 1.1
res@tmYLMinorThicknessF = 0.8
res@tmXBMinorThicknessF = 0.8
res@tmYLMajorLengthF = 0.015
res@tmXBMajorLengthF = 0.015
res@tmYLMinorLengthF = 0.011
res@tmXBMinorLengthF = 0.011
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.015
res@tmYLMajorOutwardLengthF = 0.015
res@tmXBMinorOutwardLengthF = 0.011
res@tmYLMinorOutwardLengthF = 0.011





res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 0.6
res@mpGeophysicalLineColor      = "black";"gray20"
;res@mpGridAndLimbOn = True
;res@mpGridLatSpacingF = 90
;res@mpGridLonSpacingF = 360
;res@gsnMaximize = True
res@cnFillOn             = True
res@cnLinesOn            = False
;res@cnLineThicknessF     = 6.0
;res@cnLineColor          = "red"
res@cnLineLabelsOn       = False
res@lbLabelBarOn         = False
res@lbOrientation        = "Vertical"
res@pmLabelBarWidthF     = 0.045
res@pmLabelBarHeightF    = 0.14
res@pmLabelBarOrthogonalPosF = 0.012
res@lbLabelFontHeightF  = 0.015
;res@cnMissingValFillColor = "white"
;res@gsnYRefLine = 0.0
;res@gsnYRefLineColor = "black"
;res@gsnYRefLineThicknessF = 5.0
res@cnMissingValFillColor = -1
res@cnLevelSelectionMode = "ExplicitLevels"              
; res@cnLevels             = (/140.0,155.0,170.0,185.0,200.0,215.0,230.0,245.0/)
; res@cnFillColors         = (/16,15,13,12,8,6,5,3,-1/)
res@cnFillPalette        = "precip2_17lev" 



;vcres@mpFillDrawOrder         = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.255
;vcres@vcRefAnnoSide           = "TopRight"
vcres@vcGlyphStyle            = "LineArrow"
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
vcres@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
vcres@vcMinDistanceF          = 0.03             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefLengthF            = 0.022            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.014   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 3.0            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.26 ;;正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右
vcres@vcRefAnnoString1        = "1"
vcres@vcRefMagnitudeF         = 1



sres = True
; sres@gsnAddCyclic =  True
sres@cnLineLabelsOn   = False
sres@cnConstFLabelOn  = False
sres@cnInfoLabelOn    = False 
sres@gsnDraw          = False                   
sres@gsnFrame         = False
sres@gsnLeftString    = ""
sres@gsnRightString   = ""
sres@lbLabelBarOn = False
sres@cnFillOn  = False
sres@cnLinesOn = True
sres2 = sres
sres3 = sres
; sres@cnLineColor = "red"
sres@cnLevelSelectionMode = "ExplicitLevels"
; sres@cnLevels         = (/2700.0/)
; sres@cnFillColors     = (/"red","gray20"/)
sres@cnMissingValFillColor= -1
; sres@cnLineThicknessF = 6.0
; sres@cnLineDashPattern = 14


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 4
; sres2@cnLineColors     = (/"red","red","red","red"/)
; sres2@cnLineDashPatterns = (/16,16,16,0,0,0,0/)
; sres2@cnLineThicknesses = (/4,4,4,4,4,4,4/)




cnres@tmBorderThicknessF = 2.2
cnres@gsnDraw        = False
cnres@gsnFrame       = False
cnres@lgPerimOn = False
cnres@lgLabelFontHeightF = 0.018
cnres@vpWidthF  = 0.80
cnres@vpHeightF = 0.32
cnres@vpXF      = 0.06
cnres@vpYF      = 0.5

cnres@trXMinF = -60.0
cnres@trXMaxF = 60.0
cnres@trYMinF = 0.0
cnres@trYMaxF = 85.0


cnres@tmXTOn        = False
cnres@tmYROn        = False
cnres@tmXBOn        = True
cnres@tmYLOn        = True
cnres@tmXBLabelFontHeightF = 0.0215
cnres@tmYLLabelFontHeightF = 0.0215
cnres@tmYRLabelFontHeightF = 0.0215
cnres@gsnLeftStringFontHeightF   = 0.022
cnres@gsnRightStringFontHeightF  = 0.022
cnres@gsnCenterStringFontHeightF  = 0.02
cnres@tmYRLabelsOn   = False
; cnres@tmYROn         = True


cnres@tmYLMode = "Explicit"
cnres@tmYLValues = fspan(1980.0, 2020.0, 5)
cnres@tmYLLabels = fspan(1980.0, 2020.0, 5)
cnres@tmYLMinorOn = True
cnres@tmYLMinorValues = fspan(1980.0, 2020.0, 9)

cnres@tmXBMode = "Explicit"
cnres@tmXBValues = fspan(-60, 120.0, 7)
cnres@tmXBLabels = fspan(-60, 120.0, 7)
cnres@tmXBMinorOn = True
cnres@tmXBMinorValues = fspan(-75, 120.0, 14)

cnres@cnFillOn       = True
cnres@cnLinesOn      = False
cnres@cnLineLabelsOn = False
cnres@lbLabelBarOn   = False
cnres@cnInfoLabelOn  = False
cnres@cnLevelSelectionMode = "ExplicitLevels"

cnres@tmYLMajorThicknessF = 2.5
cnres@tmXBMajorThicknessF = 2.5
cnres@tmYLMinorThicknessF = 1.5
cnres@tmXBMinorThicknessF = 1.5
cnres@tmYLMajorLengthF = 0.014
cnres@tmXBMajorLengthF = 0.014
cnres@tmYLMinorLengthF = 0.01
cnres@tmXBMinorLengthF = 0.01
cnres@tmXBLabelDeltaF = -0.4
cnres@tmYLLabelDeltaF = -0.4
cnres@tmXBMajorOutwardLengthF = 0.014
cnres@tmYLMajorOutwardLengthF = 0.014
cnres@tmXBMinorOutwardLengthF = 0.01
cnres@tmYLMinorOutwardLengthF = 0.01






lires = True
lires@tmBorderThicknessF = 1.5
lires@gsnDraw        = False
lires@gsnFrame       = False
lires@tiXAxisOn      = False
lires@tiYAxisOn      = False
lires@tmXTLabelsOn   = False
lires@tmXTOn         = False

lires@trYMinF        = 0.0
lires@trYMaxF        = 0.5
lires@trXMinF        = 0.5
lires@trXMaxF        = 12.5

lires@tmXBLabelFontHeightF = 0.042
lires@tmXTLabelFontHeightF = 0.042
lires@tmYLLabelFontHeightF = 0.042
lires@tmYRLabelFontHeightF = 0.042
lires@gsnLeftStringFontHeightF   = 0.04
lires@gsnRightStringFontHeightF  = 0.04
lires@tmYRLabelsOn   = False
lires@tmYROn         = False

lires@tmYLOn         = True
lires@tmYLLabelsOn   = True
lires@tmYLMode = "Explicit"
lires@tmXBMode = "Explicit"
; lires@tmYLValues = fspan(0, 80, 5)
; lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
; lires@tmYLMinorOn = True
; lires@tmYLMinorValues = fspan(0, 90, 10)

; lires@tmYLValues = fspan(0, 80, 5)
; lires@tmYLLabels = (/"0","20N","40N","60N","80N"/)
; lires@tmYLMinorOn = True
; lires@tmYLMinorValues = fspan(0, 90, 10)
;lires@tmYLMode = "Explicit"
;lires@tmYLValues = fspan(1, 71, 8)
;lires@tmYLLabels = (/"1-April","11-April","21-April","1-May","11-May","21-May","31-May","10-June"/)

lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
;lires@xyMarkers = (/4,1,1/)
;lires@xyMarkerSizes = (/0.015,0.001,0.001/)
;lires@xyMarkerColors= (/"black","black","black"/)
;lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)

lires@lgPerimOn = False
lires@lgLabelFontHeightF = 0.018
lires@vpWidthF  = 0.8
lires@vpHeightF = 0.7
lires@vpXF      = 0.06
lires@vpYF      = 0.75

lires@gsnLeftString  = ""
lires@gsnRightString = ""
lires@xyMarkLineModes = (/"Lines","Lines","Lines"/)
; lires@xyMarkers = (/1,1,1/)
; lires@xyMarkerSizes = (/0.015,0.001,0.001/)
; lires@xyMarkerColors= (/"black","black","black"/)
; lires@xyMarkerThicknesses = (/4.0,2.0,1.0/)
lires@xyLineThicknesses  = (/3.5,8.0,8.0/)
lires@xyLineColors       = (/"red","blue","orange"/)
lires@xyDashPatterns     = (/0,4,2/)

lires@tmYLMajorThicknessF = 2.5
lires@tmXTMajorThicknessF = 2.5
lires@tmXBMajorThicknessF = 2.5
lires@tmYLMinorThicknessF = 1.5
lires@tmXTMinorThicknessF = 1.5
lires@tmXBMinorThicknessF = 1.5
lires@tmYLMajorLengthF = 0.038
lires@tmXTMajorLengthF = 0.038
lires@tmXBMajorLengthF = 0.038
lires@tmYLMinorLengthF = 0.028
lires@tmXTMinorLengthF = 0.028
lires@tmXBMinorLengthF = 0.028
lires@tmXBLabelDeltaF = -0.25
lires@tmXTLabelDeltaF = -0.25
lires@tmYLLabelDeltaF = -0.25
lires@tmXBMajorOutwardLengthF = 0.038
lires@tmYLMajorOutwardLengthF = 0.038
lires@tmXTMajorOutwardLengthF = 0.038
lires@tmXBMinorOutwardLengthF = 0.028
lires@tmYLMinorOutwardLengthF = 0.028
lires@tmXTMinorOutwardLengthF = 0.028






;*************************打点***************************
ores                 = True   
ores@gsnAddCyclic    =  True         
ores@gsnDraw         = False               
ores@gsnFrame        = False                 
ores@cnLineLabelsOn  = False               
ores@cnLinesOn       = False                 
ores@cnInfoLabelOn   = False                                                                                                        
ores@cnFillOn        = True                               
ores@lbLabelBarOn    = False                                                                                                            
ores@cnLevelSelectionMode = "ExplicitLevels"                                                                                                                                                                                
; ores@cnLevels        = (/14.007/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ores@cnFillPatterns       = (/-1,11/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
; ores@cnFillScales         = (/0.6,0.6/)    ; change densities                                                                         
ores@cnMonoFillColor      =True                                                                                                       
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "black"

ores@cnLevels        = (/0.05/)   ;;90%
; ores@cnLevels        = (/-0.312,0.312/)   ;;95%
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
ores@cnFillPatterns       = (/-1,3/)     ; the patterns, 6                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
ores@cnFillScales         = (/0.45,0.45/)    ; change densities   
ores@cnFillDotSizeF       = 0.0025



;**********************plot**************************
plot = new(14, "graphic")
plot2 = new(20, "graphic")
plot3 = new(20, "graphic")
plot4 = new(20, "graphic")
vector = new(14, "graphic")
contour1 = new(14, "graphic")
contour2 = new(24, "graphic")
contour3 = new(14, "graphic")
contour4 = new(14, "graphic")
contour5 = new(14, "graphic")
contour6 = new(14, "graphic")
contour7 = new(14, "graphic")
contour8 = new(14, "graphic")
topooo = new(14, "graphic") 






;;-------------------------------------------------------------------------------------------
pltType = "eps"
pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure1a"
wks = gsn_open_wks(pltType, pltName)





res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/14,17,20,23,26,29/)
res@cnFillColors := (/-1,42,84,126,168,205,225/)
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.3
res@mpGeophysicalLineColor = "black"
plot(0) = gsn_csm_contour_map_ce(wks, wbt95_clm, res)
plot2(0) = gsn_csm_contour(wks, both_high, ores)
overlay(plot(0), plot2(0))
; plot2(1) = gsn_csm_contour(wks, both_high, ores)
; overlay(plot(1), plot2(1))



res@cnFillPalette = "WhiteBlueGreenYellowRed" ; "WhiteBlue"
res@cnLevels     := (/17,22,27,32,37,42/)
res@cnFillColors := (/-1,42,84,126,168,205,225/)
res@tmXBLabelsOn  = True
res@tmYLLabelsOn  = True
; res@mpGeophysicalLineThicknessF = 0.8
plot(1) = gsn_csm_contour_map_ce(wks, degree, res)
plot2(1) = gsn_csm_contour(wks, both_high, ores)
overlay(plot(1), plot2(1))




res@tmYLLabelsOn  = False
res@cnLevels     := (/17,22,27,32,37,42/)
plot(2) = gsn_csm_contour_map_ce(wks, degree_HWs_Wet, res)

res@tmYLLabelsOn  = False
; res@cnLevels     := (/15,21,27,33,39,45/)
; plot(3) = gsn_csm_contour_map_ce(wks, degree_Wet_Wet, res)

res@tmYLLabelsOn  = True
res@cnFillPalette = "WhiteGreen"
res@cnLevels     := (/0.6,0.65,0.7,0.75,0.8,0.85/)
res@cnFillColors := (/-1,40,70,100,130,160,190/)
plot(3) = gsn_csm_contour_map_ce(wks, Cor, res)




  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation        = "Vertical"
  pres@pmLabelBarWidthF     = 0.035
  pres@pmLabelBarHeightF    = 0.15
  pres@pmLabelBarOrthogonalPosF = 0.003
  ; pres@pmLabelBarParallelPosF = 0.01
  pres@lbLabelFontHeightF  = 0.008

  pres@gsnPanelTop  = 0.71
  pres@gsnPanelBottom = 0.51
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.525
  ; pres@amJust  ="TopLeft"
  ; pres@gsnPanelFigureStringsFontHeightF = 0.016
;   pres@gsnPanelFigureStringsPerimOn = False
; ;   pres@gsnPanelFigureStringsBackgroundFillColor = -1
;   pres@gsnPanelFigureStrings = (/"","(c)"/)
  gsn_panel(wks,plot(0),(/1,1/),pres)

  pres@gsnPanelTop  = 0.467
  pres@gsnPanelBottom = 0.267
  pres@gsnPanelLeft  = 0.01
  pres@gsnPanelRight = 0.525
  gsn_panel(wks,plot(1),(/1,1/),pres)


  pres@gsnPanelTop  = 0.467-0.010
  pres@gsnPanelBottom = 0.267+0.006
  pres@gsnPanelLeft  = 0.495
  pres@gsnPanelRight = 0.99
  gsn_panel(wks,plot(2),(/1,1/),pres)  

  ; pres@gsnPanelTop  = 0.465
  ; pres@gsnPanelBottom = 0.27
  ; pres@gsnPanelLeft  = 0.495
  ; pres@gsnPanelRight = 0.99
  ; gsn_panel(wks,plot(3),(/1,1/),pres)




txres               = True
txres@txFontColor = "black"
txres@txFontHeightF = 0.0105
txres@txFont = "helvetica"
gsn_text_ndc(wks, "Wet   bulb   temperature",  0.235, 0.715, txres)
gsn_text_ndc(wks, "Proportion   of   significant   links",  0.715, 0.715, txres)
gsn_text_ndc(wks, "Spatial   aggregation   of   humid   heatwaves",  0.235, 0.475, txres)
gsn_text_ndc(wks, "Spatial   concurrence   with   wet   events",  0.697, 0.475, txres)

; gsn_text_ndc(wks, "Local   concurrences   of   humid   heatwaves",  0.697, 0.715, txres)
; gsn_text_ndc(wks, "Local   concurrences   with   wet   events",  0.235, 0.475, txres)
; gsn_text_ndc(wks, "Proportion   of   significant   links",  0.715, 0.475, txres)

txres@txFontHeightF = 0.0135
txres@txFont = "helvetica-bold"
gsn_text_ndc(wks, "a",  0.03, 0.73, txres)
gsn_text_ndc(wks, "b",  0.52, 0.73, txres)
gsn_text_ndc(wks, "c",  0.03, 0.49, txres)
gsn_text_ndc(wks, "d",  0.51, 0.49, txres)


frame(wks)
   delete(wks)  
   system("convert -geometry 2400x2400 -density 1600x1600 -trim " + pltName + "."+pltType + " " + pltName + ".png")     
  system("rm " + pltName + "." + pltType + " -f")

end







