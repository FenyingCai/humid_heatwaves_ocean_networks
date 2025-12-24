load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
;;----------------------------------------------;;
;;--------------    read network   -------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part4_wbt95_qs95_SH/wbt95_qs95_SH2x2_1982_2023_sign99.nc","r")
networks0 = f1->networks0(:,:,:,:)
networks0(0:14,:,:,:) = 1.0e20  ;; exclude 90S-60S
networks0(:,:,0:14,:) = 1.0e20
networks0@_FillValue = 1.0e20
dim_0 = dimsizes(networks0)


;;--------------- remove ocean -----------------;;
f00 = addfile("/public/home/fcai/data0/topography/topo_2x2.nc","r")
topo = f00->topo({0:90},:)
topo_networks1 = conform(networks0, topo, (/0,1/))
topo_networks2 = conform(networks0, topo, (/2,3/))

; networks0 = where(topo_networks1.ge.0.0, networks0, 1.0e20)
; networks0 = where(topo_networks2.lt.0.0, networks0, 1.0e20)
networks0@_FillValue = 1.0e20
printMinMax(networks0, 0)
networksRef = where(topo_networks1.ge.-9999.0, 1.0, 1.0e20)
; networksRef = where(topo_networks2.lt.0.0, networksRef, 1.0e20)
networksRef@_FillValue = 1.0e20
print("all link number  =  "+dim_sum(dim_sum(dim_sum(dim_sum(networksRef)))))



networks0 = where(networks0.gt.0.0, 1.0, 1.0e20)
networks0@_FillValue = 1.0e20
print("significant link number  =  "+dim_sum(dim_sum(dim_sum(dim_sum(networks0)))))
; networks0 = where(networks0.ge.41.0, 1.0, 1.0e20)
; networks0@_FillValue = 1.0e20
; print("significant link number (>41)  =  "+dim_sum(dim_sum(dim_sum(dim_sum(networks0)))))





;;--------------------------  lat weighted  ---------------------------;;
lat_matrix = new((/dim_0(0),dim_0(1), dim_0(0),dim_0(1)/), float)
do i = 0,dim_0(0)-1
  do j = 0,dim_0(0)-1
    lat_matrix(i,:,j,:) = cos((i*2-89+j*2-89) /2.0 /180.0 * 3.14159)
  end do
end do




;;----------------------------------------------  read distance  ---------------------------------------------;;
f0 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/1_threshold_2x2/distance_i_j_2x2_SH.nc","r")
distance = f0->distance(:,:,:,:)
distance(0:14,:,:,:) = 1.0e20  ;; exclude 90S-60S
distance(:,:,0:14,:) = 1.0e20
distance@_FillValue = 1.0e20
dist = fspan(0.0,10200.0,52)


pdf = new(52, float)
ref = new(52, float)
do i =0,51

if (dist(i).eq.10200.0) then
networks0_new = where((distance.ge.dist(i)-100), networks0, 0.0)
networksRef_new = where((distance.ge.dist(i)-100), networksRef, 0.0)
pdf(i) = dim_sum(dim_sum(dim_sum(dim_sum(networks0_new*lat_matrix))))
ref(i) = dim_sum(dim_sum(dim_sum(dim_sum(networksRef_new*lat_matrix))))
end if
if (dist(i).lt.10200.0) then
networks0_new = where((distance.ge.dist(i)-100 .and. distance.lt.dist(i)+100), networks0, 0.0)
networksRef_new = where((distance.ge.dist(i)-100 .and. distance.lt.dist(i)+100), networksRef, 0.0)
pdf(i) = dim_sum(dim_sum(dim_sum(dim_sum(networks0_new*lat_matrix))))
ref(i) = dim_sum(dim_sum(dim_sum(dim_sum(networksRef_new*lat_matrix))))
end if

delete(networks0_new)
print(" i="+i+",  dist= "+dist(i)+",  prob = "+pdf(i)/ref(i))
end do

delete([/lat_matrix,distance,networks0/])
print(pdf/ref)

; pdf = pdf/dim_sum(pdf)
pdf_ref = pdf/ref


print(dist)
print("----  Land-Ocean,  concurrent days > 99% significant  ----")
print(pdf_ref)




;;---------------------------------- to nc file -----------------------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/2_network/part4_wbt95_qs95_SH/"
name1="distance_PDF_1982_2023_wbt95_qs95_SH"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->pdf=pdf
ncdf1->ref=ref
ncdf1->pdf_ref=pdf_ref


end






