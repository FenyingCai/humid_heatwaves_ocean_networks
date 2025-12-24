load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin


;;------------------------------------------------------------------------------------------------------;;
;;------------------------------------------------  NH  ------------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 := (f1->topo({-90:90},:))
dim_0 = dimsizes(topo0)
lat = topo0&lat
lon = topo0&lon
printVarSummary(topo0)



;;------------------------------------------- read composite -------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/6_response/composite_analysis/SAsia_composite.nc","r")
q_SAsia_diff=f1->q_diff
u_SAsia_diff=f1->u_diff
v_SAsia_diff=f1->v_diff
qu_SAsia_diff=f1->qu_diff
qv_SAsia_diff=f1->qv_diff
sst_SAsia_diff=f1->sst_diff
wbt_SAsia_diff=f1->wbt_diff

q_SAsia_t=f1->q_t
u_SAsia_t=f1->u_t
v_SAsia_t=f1->v_t
qu_SAsia_t=f1->qu_t
qv_SAsia_t=f1->qv_t
sst_SAsia_t=f1->sst_t
wbt_SAsia_t=f1->wbt_t

qu_SAsia_diff = qu_SAsia_diff*1000.0
qv_SAsia_diff = qv_SAsia_diff*1000.0
qu_SAsia_diff = where(qu_SAsia_t.le.0.05.or.qv_SAsia_t.le.0.05, qu_SAsia_diff, 1.0e20)
qv_SAsia_diff = where(qu_SAsia_t.le.0.05.or.qv_SAsia_t.le.0.05, qv_SAsia_diff, 1.0e20)
qu_SAsia_diff@_FillValue = 1.0e20
qv_SAsia_diff@_FillValue = 1.0e20
copy_VarCoords(u_SAsia_diff, qu_SAsia_diff)
copy_VarCoords(v_SAsia_diff, qv_SAsia_diff)
q_SAsia_diff = q_SAsia_diff*1000.0
copy_VarCoords(u_SAsia_diff, q_SAsia_diff)

; q_SAsia_diff = where(q_SAsia_t.le.0.05, q_SAsia_diff, 1.0e20)
q_SAsia_diff@_FillValue = 1.0e20
copy_VarCoords(q_SAsia_t, q_SAsia_diff)
; sst_SAsia_diff = where(sst_SAsia_t.le.0.05, sst_SAsia_diff, 1.0e20)
; wbt_SAsia_diff = where(wbt_SAsia_t.le.0.05, wbt_SAsia_diff, 1.0e20)
sst_SAsia_diff@_FillValue = 1.0e20
wbt_SAsia_diff@_FillValue = 1.0e20
copy_VarCoords(sst_SAsia_t, sst_SAsia_diff)
copy_VarCoords(sst_SAsia_t, wbt_SAsia_diff)









;;------------------------------------------- read CESM output -------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/6_response/cesm/cesm_CTL_JJA.nc","r")
TSMX_ctl = f1->TSMX
QREFHT_ctl = f1->QREFHT
TREFHT_ctl = f1->TREFHT
PS_ctl = f1->PS0
U_ctl = f1->U(:,{850},:,:)
V_ctl = f1->V(:,{850},:,:)
Q_ctl = f1->Q(:,{850},:,:)
Q_ctl = Q_ctl * 1000.0
QU_ctl = Q_ctl * U_ctl
QV_ctl = Q_ctl * V_ctl
copy_VarCoords(Q_ctl, QU_ctl)
copy_VarCoords(Q_ctl, QV_ctl)
printVarSummary(TSMX_ctl)


f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/6_response/cesm/cesm_South_Asia2_JJA.nc","r")
TSMX_NIO = f1->TSMX
QREFHT_NIO = f1->QREFHT
TREFHT_NIO = f1->TREFHT
PS_NIO = f1->PS0
U_NIO = f1->U(:,{850},:,:)
V_NIO = f1->V(:,{850},:,:)
Q_NIO = f1->Q(:,{850},:,:)
Q_NIO = Q_NIO * 1000.0
QU_NIO = Q_NIO * U_NIO
QV_NIO = Q_NIO * V_NIO
copy_VarCoords(Q_ctl, QU_NIO)
copy_VarCoords(Q_ctl, QV_NIO)


;;--------------- calculate wet-bulb temperature ----------------;;
print(" -------------- ps --------------- ")
printMinMax(PS_ctl, 1)  ;; Pa
print(" -------------- q --------------- ")
printMinMax(QREFHT_ctl, 1)  ;; kg/kg
print(" -------------- temp --------------- ")
printMinMax(TREFHT_ctl, 1)  ;; K

RH_ctl = relhum(TREFHT_ctl, QREFHT_ctl, PS_ctl)  ;; Pa, K, kg/kg -> %
dewT_ctl = dewtemp_trh(TREFHT_ctl, RH_ctl)  ;; K, % -> K
PS_ctl = PS_ctl / 1000.0
TREFHT_ctl = TREFHT_ctl - 273.15
dewT_ctl = dewT_ctl - 273.15
print(" -------------- second --------------- ")
printMinMax(PS_ctl, 1) 
printMinMax(TREFHT_ctl, 1) 
printMinMax(dewT_ctl, 1) 
wbt_ctl = wetbulb(PS_ctl, TREFHT_ctl, dewT_ctl) ;; hPa, C, C -> C
copy_VarCoords(PS_ctl, wbt_ctl)


;;--------------- calculate wet-bulb temperature ----------------;;
print(" -------------- ps --------------- ")
printMinMax(PS_NIO, 1)  ;; Pa
print(" -------------- q --------------- ")
printMinMax(QREFHT_NIO, 1)  ;; kg/kg
print(" -------------- temp --------------- ")
printMinMax(TREFHT_NIO, 1)  ;; K

RH_NIO = relhum(TREFHT_NIO, QREFHT_NIO, PS_NIO)  ;; Pa, K, kg/kg -> %
dewT_NIO = dewtemp_trh(TREFHT_NIO, RH_NIO)  ;; K, % -> K
PS_NIO = PS_NIO / 1000.0
TREFHT_NIO = TREFHT_NIO - 273.15
dewT_NIO = dewT_NIO - 273.15
print(" -------------- second --------------- ")
printMinMax(PS_NIO, 1) 
printMinMax(TREFHT_NIO, 1) 
printMinMax(dewT_NIO, 1) 
wbt_NIO = wetbulb(PS_NIO, TREFHT_NIO, dewT_NIO) ;; hPa, C, C -> C
copy_VarCoords(PS_NIO, wbt_NIO)


;;--------------- calculate difference ----------------;;
QREFHT_ctl = QREFHT_ctl * 1000.0
QREFHT_NIO = QREFHT_NIO * 1000.0
TSMX_NIO_diff = TSMX_NIO - TSMX_ctl
QREFHT_NIO_diff = QREFHT_NIO - QREFHT_ctl
TREFHT_NIO_diff = TREFHT_NIO - TREFHT_ctl
wbt_NIO_diff = wbt_NIO - wbt_ctl
U_NIO_diff = U_NIO - U_ctl
V_NIO_diff = V_NIO - V_ctl
Q_NIO_diff = Q_NIO - Q_ctl
QU_NIO_diff = QU_NIO - QU_ctl
QV_NIO_diff = QV_NIO - QV_ctl
copy_VarCoords(TSMX_ctl, TSMX_NIO_diff)
copy_VarCoords(TSMX_ctl, QREFHT_NIO_diff)
copy_VarCoords(TSMX_ctl, TREFHT_NIO_diff)
copy_VarCoords(TSMX_ctl, wbt_NIO_diff)
copy_VarCoords(TSMX_ctl, U_NIO_diff)
copy_VarCoords(TSMX_ctl, V_NIO_diff)
copy_VarCoords(TSMX_ctl, Q_NIO_diff)
copy_VarCoords(TSMX_ctl, QU_NIO_diff)
copy_VarCoords(TSMX_ctl, QV_NIO_diff)


TSMX_NIO_diff_avg = dim_avg_n_Wrap(TSMX_NIO_diff, 0)
QREFHT_NIO_diff_avg = dim_avg_n_Wrap(QREFHT_NIO_diff, 0)
TREFHT_NIO_diff_avg = dim_avg_n_Wrap(TREFHT_NIO_diff, 0)
wbt_NIO_diff_avg = dim_avg_n_Wrap(wbt_NIO_diff, 0)
U_NIO_diff_avg = dim_avg_n_Wrap(U_NIO_diff, 0)
V_NIO_diff_avg = dim_avg_n_Wrap(V_NIO_diff, 0)
Q_NIO_diff_avg = dim_avg_n_Wrap(Q_NIO_diff, 0)
QU_NIO_diff_avg = dim_avg_n_Wrap(QU_NIO_diff, 0)
QV_NIO_diff_avg = dim_avg_n_Wrap(QV_NIO_diff, 0)



;;-------------------------- remove topo (land area & Tibetan Plateau) ----------------------------;;
topo0_cesm = linint2(topo0&lon,topo0&lat,topo0, True,TSMX_NIO_diff_avg&lon,TSMX_NIO_diff_avg&lat, 0)
TSMX_NIO_diff_avg = where(topo0_cesm.le.-80.0, TSMX_NIO_diff_avg, 0.0)
copy_VarCoords(QREFHT_NIO_diff_avg, TSMX_NIO_diff_avg)
Q_NIO_diff_avg = where(topo0_cesm.le.2000.0, Q_NIO_diff_avg, 1.0e20)
QU_NIO_diff_avg = where(topo0_cesm.le.2000.0, QU_NIO_diff_avg, 1.0e20)
QV_NIO_diff_avg = where(topo0_cesm.le.2000.0, QV_NIO_diff_avg, 1.0e20)
Q_NIO_diff_avg@_FillValue = 1.0e20
QU_NIO_diff_avg@_FillValue = 1.0e20
QV_NIO_diff_avg@_FillValue = 1.0e20
copy_VarCoords(QREFHT_NIO_diff_avg, Q_NIO_diff_avg)
copy_VarCoords(QREFHT_NIO_diff_avg, QU_NIO_diff_avg)
copy_VarCoords(QREFHT_NIO_diff_avg, QV_NIO_diff_avg)



;;------------------------- p-values ----------------------------;;
TSMX_ctl_clm = dim_avg_n_Wrap(TSMX_ctl, 0)
TSMX_ctl_var = dim_variance_n_Wrap(TSMX_ctl, 0)
TSMX_ctl_size = equiv_sample_size(TSMX_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(TSMX_ctl(0,:,:), TSMX_ctl_size)
QREFHT_ctl_clm = dim_avg_n_Wrap(QREFHT_ctl, 0)
QREFHT_ctl_var = dim_variance_n_Wrap(QREFHT_ctl, 0)
QREFHT_ctl_size = equiv_sample_size(QREFHT_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(QREFHT_ctl(0,:,:), QREFHT_ctl_size)
TREFHT_ctl_clm = dim_avg_n_Wrap(TREFHT_ctl, 0)
TREFHT_ctl_var = dim_variance_n_Wrap(TREFHT_ctl, 0)
TREFHT_ctl_size = equiv_sample_size(TREFHT_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(TREFHT_ctl(0,:,:), TREFHT_ctl_size)
wbt_ctl_clm = dim_avg_n_Wrap(wbt_ctl, 0)
wbt_ctl_var = dim_variance_n_Wrap(wbt_ctl, 0)
wbt_ctl!0="year"
wbt_ctl_size = equiv_sample_size(wbt_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(wbt_ctl(0,:,:), wbt_ctl_size)
U_ctl_clm = dim_avg_n_Wrap(U_ctl, 0)
U_ctl_var = dim_variance_n_Wrap(U_ctl, 0)
U_ctl_size = equiv_sample_size(U_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(U_ctl(0,:,:), U_ctl_size)
V_ctl_clm = dim_avg_n_Wrap(V_ctl, 0)
V_ctl_var = dim_variance_n_Wrap(V_ctl, 0)
V_ctl_size = equiv_sample_size(V_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(V_ctl(0,:,:), V_ctl_size)
Q_ctl_clm = dim_avg_n_Wrap(Q_ctl, 0)
Q_ctl_var = dim_variance_n_Wrap(Q_ctl, 0)
Q_ctl_size = equiv_sample_size(Q_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(Q_ctl(0,:,:), Q_ctl_size)
QU_ctl_clm = dim_avg_n_Wrap(QU_ctl, 0)
QU_ctl_var = dim_variance_n_Wrap(QU_ctl, 0)
QU_ctl_size = equiv_sample_size(QU_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(QU_ctl(0,:,:), QU_ctl_size)
QV_ctl_clm = dim_avg_n_Wrap(QV_ctl, 0)
QV_ctl_var = dim_variance_n_Wrap(QV_ctl, 0)
QV_ctl_size = equiv_sample_size(QV_ctl(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(QV_ctl(0,:,:), QV_ctl_size)

TSMX_NIO_clm = dim_avg_n_Wrap(TSMX_NIO, 0)
TSMX_NIO_var = dim_variance_n_Wrap(TSMX_NIO, 0)
TSMX_NIO_size = equiv_sample_size(TSMX_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(TSMX_NIO(0,:,:), TSMX_NIO_size)
QREFHT_NIO_clm = dim_avg_n_Wrap(QREFHT_NIO, 0)
QREFHT_NIO_var = dim_variance_n_Wrap(QREFHT_NIO, 0)
QREFHT_NIO_size = equiv_sample_size(QREFHT_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(QREFHT_NIO(0,:,:), QREFHT_NIO_size)
TREFHT_NIO_clm = dim_avg_n_Wrap(TREFHT_NIO, 0)
TREFHT_NIO_var = dim_variance_n_Wrap(TREFHT_NIO, 0)
TREFHT_NIO_size = equiv_sample_size(TREFHT_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(TREFHT_NIO(0,:,:), TREFHT_NIO_size)
wbt_NIO_clm = dim_avg_n_Wrap(wbt_NIO, 0)
wbt_NIO_var = dim_variance_n_Wrap(wbt_NIO, 0)
wbt_NIO!0="year"
wbt_NIO_size = equiv_sample_size(wbt_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(wbt_NIO(0,:,:), wbt_NIO_size)
U_NIO_clm = dim_avg_n_Wrap(U_NIO, 0)
U_NIO_var = dim_variance_n_Wrap(U_NIO, 0)
U_NIO_size = equiv_sample_size(U_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(U_NIO(0,:,:), U_NIO_size)
V_NIO_clm = dim_avg_n_Wrap(V_NIO, 0)
V_NIO_var = dim_variance_n_Wrap(V_NIO, 0)
V_NIO_size = equiv_sample_size(V_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(V_NIO(0,:,:), V_NIO_size)
Q_NIO_clm = dim_avg_n_Wrap(Q_NIO, 0)
Q_NIO_var = dim_variance_n_Wrap(Q_NIO, 0)
Q_NIO_size = equiv_sample_size(Q_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(Q_NIO(0,:,:), Q_NIO_size)
QU_NIO_clm = dim_avg_n_Wrap(QU_NIO, 0)
QU_NIO_var = dim_variance_n_Wrap(QU_NIO, 0)
QU_NIO_size = equiv_sample_size(QU_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(QU_NIO(0,:,:), QU_NIO_size)
QV_NIO_clm = dim_avg_n_Wrap(QV_NIO, 0)
QV_NIO_var = dim_variance_n_Wrap(QV_NIO, 0)
QV_NIO_size = equiv_sample_size(QV_NIO(lat|:,lon|:,year|:), 0.05, 0)
copy_VarCoords(QV_NIO(0,:,:), QV_NIO_size)


TSMX_NIO_t = ttest(TSMX_NIO_clm, TSMX_NIO_var, TSMX_NIO_size, TSMX_ctl_clm, TSMX_ctl_var, TSMX_ctl_size, False, False)
copy_VarCoords(TSMX_NIO_clm, TSMX_NIO_t)
QREFHT_NIO_t = ttest(QREFHT_NIO_clm, QREFHT_NIO_var, QREFHT_NIO_size, QREFHT_ctl_clm, QREFHT_ctl_var, QREFHT_ctl_size, False, False)
copy_VarCoords(QREFHT_NIO_clm, QREFHT_NIO_t)
TREFHT_NIO_t = ttest(TREFHT_NIO_clm, TREFHT_NIO_var, TREFHT_NIO_size, TREFHT_ctl_clm, TREFHT_ctl_var, TREFHT_ctl_size, False, False)
copy_VarCoords(TREFHT_NIO_clm, TREFHT_NIO_t)
wbt_NIO_t = ttest(wbt_NIO_clm, wbt_NIO_var, wbt_NIO_size, wbt_ctl_clm, wbt_ctl_var, wbt_ctl_size, False, False)
copy_VarCoords(wbt_NIO_clm, wbt_NIO_t)
U_NIO_t = ttest(U_NIO_clm, U_NIO_var, U_NIO_size, U_ctl_clm, U_ctl_var, U_ctl_size, False, False)
copy_VarCoords(U_NIO_clm, U_NIO_t)
V_NIO_t = ttest(V_NIO_clm, V_NIO_var, V_NIO_size, V_ctl_clm, V_ctl_var, V_ctl_size, False, False)
copy_VarCoords(V_NIO_clm, V_NIO_t)
Q_NIO_t = ttest(Q_NIO_clm, Q_NIO_var, Q_NIO_size, Q_ctl_clm, Q_ctl_var, Q_ctl_size, False, False)
copy_VarCoords(Q_NIO_clm, Q_NIO_t)
QU_NIO_t = ttest(QU_NIO_clm, QU_NIO_var, QU_NIO_size, QU_ctl_clm, QU_ctl_var, QU_ctl_size, False, False)
copy_VarCoords(QU_NIO_clm, QU_NIO_t)
QV_NIO_t = ttest(QV_NIO_clm, QV_NIO_var, QV_NIO_size, QV_ctl_clm, QV_ctl_var, QV_ctl_size, False, False)
copy_VarCoords(QV_NIO_clm, QV_NIO_t)

QU_NIO_diff_avg = where(QU_NIO_t.le.0.05.or.QV_NIO_t.le.0.05, QU_NIO_diff_avg, 1.0e20)
QV_NIO_diff_avg = where(QU_NIO_t.le.0.05.or.QV_NIO_t.le.0.05, QV_NIO_diff_avg, 1.0e20)
QU_NIO_diff_avg@_FillValue = 1.0e20
QV_NIO_diff_avg@_FillValue = 1.0e20
copy_VarCoords(QREFHT_NIO_diff_avg, QU_NIO_diff_avg)
copy_VarCoords(QREFHT_NIO_diff_avg, QV_NIO_diff_avg)

wbt_NIO_diff_avg = where(topo0_cesm.gt.0.0, wbt_NIO_diff_avg, 1.0e20)
wbt_NIO_t = where(topo0_cesm.gt.0.0, wbt_NIO_t, 1.0e20)
wbt_NIO_diff_avg@_FillValue = 1.0e20
wbt_NIO_t@_FillValue = 1.0e20
copy_VarCoords(QREFHT_NIO_diff_avg, wbt_NIO_diff_avg)
copy_VarCoords(QREFHT_NIO_diff_avg, wbt_NIO_t)












;;------------------------------  read lat, lon  ------------------------------;;
f1 = addfile("/public/data2/model/fcai/data0/topography/topo_adapt_era5_air2m.nc", "r")
topo0 := (f1->topo({0:90},{0:180}))
topo_NHSH := (f1->topo({-30:60},{0:180}))
lat := topo0&lat
lon := topo0&lon
lat_NHSH := topo_NHSH&lat

lat_weight = new((/30*33, 91,181/), float)
do ilat = 0,90
lat_weight(:,ilat,:) = cos(abs(lat_NHSH(ilat)) / 180.0 * 3.14159)
end do
printMinMax(lat_weight, 1)

topo3D = conform(lat_weight, topo_NHSH, (/1,2/))
lat_weight_land = where(topo3D.ge.0.0, lat_weight, 0.0)
lat_weight_ocean = where(topo3D.lt.0.0, lat_weight, 0.0)

lat_weight_land!1 = "lat"
lat_weight_land!2 = "lon"
lat_weight_land&lat = lat_NHSH
lat_weight_land&lon = lon
copy_VarCoords(lat_weight_land, lat_weight_ocean)


;;---------------------------------------------------------------------------------------------------;;
f1 = addfile("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig5_CMIP6_historical/wbt_historical_30models_1979_2014_30S60N.nc","r")
wbt0 = dim_avg_n_Wrap(f1->wbt(:,3:35,:,:,0:180), 2)  ;; 1982-2023
wbt0 = dtrend_msg_n(fspan(1, 33, 33), wbt0, True, True, 1)
dim_wbt = dimsizes(wbt0)
wbt = new((/33*30,dim_wbt(2),dim_wbt(3)/), float)
do imodel = 0,29
wbt(imodel*33:imodel*33+32,:,:) = wbt0(imodel,:,:,:)
end do


topo_3D = conform(wbt, topo_NHSH, (/1,2/))
wbt = where(topo_3D.ge.-50.0, wbt, 1.0e20)
wbt@_FillValue = 1.0e20
wbt!1 = "lat"
wbt!2 = "lon"
wbt&lat = lat_NHSH
wbt&lon = lon
printVarSummary(wbt)

SA_lat1 = 8
SA_lat2 = 28
SA_lon1 = 70
SA_lon2 = 90
wbt = where(wbt.ge.-100.0, wbt, 0.0)
SA_wbt = dim_sum_n_Wrap(dim_sum_n_Wrap(wbt(:,{SA_lat1:SA_lat2},{SA_lon1:SA_lon2})*lat_weight_land(:,{SA_lat1:SA_lat2},{SA_lon1:SA_lon2}), 2), 1) / dim_sum_n_Wrap(dim_sum_n_Wrap(lat_weight_land(:,{SA_lat1:SA_lat2},{SA_lon1:SA_lon2}), 2), 1)
; print(SA_wbt)

;;------------------------------------------------------------------------------------------------------------------------;;
f1 = addfile("/public/data2/model/fcai/CMIP6/historical_tos_mon/tos_for8humidheat/historical_tos_12mons_30models_1979_2023.nc","r")
sst0 = doubletofloat(f1->tos(:,3:35,5:7,60:150,0:180))  ;; 1982-2023
sst0@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(sst0))) then
      replace_ieeenan (sst0, value, 0)
      sst0@_FillValue = value
end if
sst0_ = dim_avg_n_Wrap(sst0, 2)
sst0_ = dtrend_msg_n(fspan(1, 33, 33), sst0_, True, True, 1)
dim_sst = dimsizes(sst0)
sst = new((/33*30,dim_sst(3),dim_sst(4)/), float)
do imodel = 0,29
sst(imodel*33:imodel*33+32,:,:) = sst0_(imodel,:,:,:)
end do
sst!1 = "lat"
sst!2 = "lon"
sst&lat = topo_NHSH&lat
sst&lon = lon
NIO_lat1 = 0
NIO_lat2 = 25
NIO_lon1 = 50
NIO_lon2 = 100
NIO_sst = dim_sum_n_Wrap(dim_sum_n_Wrap(sst(:,{NIO_lat1:NIO_lat2},{NIO_lon1:NIO_lon2})*lat_weight_ocean(:,{NIO_lat1:NIO_lat2},{NIO_lon1:NIO_lon2}), 2), 1) / dim_sum_n_Wrap(dim_sum_n_Wrap(lat_weight_ocean(:,{NIO_lat1:NIO_lat2},{NIO_lon1:NIO_lon2}), 2), 1)


f2 = addfile("/public/data2/model/fcai/CMIP6/historical_hus_mon/hus_for8humidheat/historical_hus850_JJA_30models_1979_2014.nc","r")
q0 = doubletofloat(f2->hus(:,3:35,60:150,0:180))  ;; 1982-2023
q0@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(q0))) then
      replace_ieeenan (q0, value, 0)
      q0@_FillValue = value
end if
q0 = dtrend_msg_n(fspan(1, 33, 33), q0, False, False, 1)
q0@_FillValue = 1.0e20

f2 = addfile("/public/data2/model/fcai/CMIP6/historical_ua_mon/ua_for8humidheat/historical_ua850_JJA_30models_1979_2014.nc","r")
u0 = doubletofloat(f2->ua(:,3:35,60:150,0:180))  ;; 1982-2023
u0@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(u0))) then
      replace_ieeenan (u0, value, 0)
      u0@_FillValue = value
end if
u0 = dtrend_msg_n(fspan(1, 33, 33), u0, False, False, 1)
q0@_FillValue = 1.0e20

f3 = addfile("/public/data2/model/fcai/CMIP6/historical_va_mon/va_for8humidheat/historical_va850_JJA_30models_1979_2014.nc","r")
v0 = doubletofloat(f3->va(:,3:35,60:150,0:180))  ;; 1982-2023
v0@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(v0))) then
      replace_ieeenan (v0, value, 0)
      v0@_FillValue = value
end if
v0 = dtrend_msg_n(fspan(1, 33, 33), v0, False, False, 1)
q0@_FillValue = 1.0e20


q = new((/33*30,dim_sst(3),dim_sst(4)/), float)
qu = new((/33*30,dim_sst(3),dim_sst(4)/), float)
qv = new((/33*30,dim_sst(3),dim_sst(4)/), float)
do imodel = 0,29
q(imodel*33:imodel*33+32,:,:) = q0(imodel,:,:,:) * 1000.0
qu(imodel*33:imodel*33+32,:,:) = q0(imodel,:,:,:) * u0(imodel,:,:,:) * 1000.0
qv(imodel*33:imodel*33+32,:,:) = q0(imodel,:,:,:) * v0(imodel,:,:,:) * 1000.0
end do

SA_wbt_= dim_standardize_n_Wrap(SA_wbt, 1, 0)
SA_corq = escorc_n(SA_wbt_, q, 0, 0)
SA_corqu = escorc_n(SA_wbt_, qu, 0, 0)
SA_corqv = escorc_n(SA_wbt_, qv, 0, 0)
SA_corwbt = escorc_n(SA_wbt_, wbt, 0, 0)
SA_regq = regCoef_n(SA_wbt_, q, 0, 0) * 100.0
SA_regqu = regCoef_n(SA_wbt_, qu, 0, 0)
SA_regqv = regCoef_n(SA_wbt_, qv, 0, 0)
SA_regwbt = regCoef_n(SA_wbt_, wbt, 0, 0) * 100.0

SA_corSST = escorc_n(SA_wbt, sst, 0, 0)
copy_VarCoords(topo_NHSH, SA_corSST)
SA_regSST = regCoef_n(SA_wbt_, sst, 0, 0) * 100.0
copy_VarCoords(topo_NHSH, SA_regSST)

; SA_regq = where(abs(SA_corq).ge.0.101, SA_regq, 1.0e20)
SA_regq@_FillValue = 1.0e20
SA_regqu = where(abs(SA_corqu).ge.0.101.or.abs(SA_corqv).ge.0.101, SA_regqu, 1.0e20)
SA_regqv = where(abs(SA_corqu).ge.0.101.or.abs(SA_corqv).ge.0.101, SA_regqv, 1.0e20)
SA_regq = where(topo_NHSH.le.2000.0, SA_regq, 1.0e20)
SA_regqu = where(topo_NHSH.le.2000.0, SA_regqu, 1.0e20)
SA_regqv = where(topo_NHSH.le.2000.0, SA_regqv, 1.0e20)
SA_regq@_FillValue = 1.0e20
SA_regqu@_FillValue = 1.0e20
SA_regqv@_FillValue = 1.0e20
SA_regwbt = where(topo_NHSH.gt.0.0, SA_regwbt, 1.0e20)
SA_corwbt = where(topo_NHSH.gt.0.0, SA_corwbt, 1.0e20)
SA_corwbt = where(abs(SA_regwbt).ge.10.0, SA_corwbt, 1.0e20)
SA_regwbt@_FillValue = 1.0e20
SA_corwbt@_FillValue = 1.0e20
copy_VarCoords(topo_NHSH, SA_corq)
copy_VarCoords(topo_NHSH, SA_corqu)
copy_VarCoords(topo_NHSH, SA_corqv)
copy_VarCoords(wbt(0,:,:), SA_corwbt)
copy_VarCoords(topo_NHSH, SA_regq)
copy_VarCoords(topo_NHSH, SA_regqu)
copy_VarCoords(topo_NHSH, SA_regqv)
copy_VarCoords(wbt(0,:,:), SA_regwbt)


;;---------------------------------------------------;;
key_region3 = topo0
key_region3 = 0.0
key_region3({SA_lat1:SA_lat2},{SA_lon1:SA_lon2}) = 1.0
key_region3 = where(topo0.ge.0.0, key_region3, -1.0)
copy_VarCoords(topo0, key_region3)






;;************************************************** plots **********************************************;;
res = True
res@tmBorderThicknessF = 1.2
res@mpPerimOn = False
; res@tmXTBorderOn = False
; res@tmXBBorderOn = False
; res@tmYLBorderOn = False
; res@tmYRBorderOn = False
res@gsnDraw      =  False
res@gsnFrame     =  False
res@gsnAddCyclic =  False
res@gsnRightString       = ""
res@gsnLeftString        = ""
res@gsnLeftStringFontHeightF   = 0.04
res@gsnRightStringFontHeightF  = 0.04
res@gsnCenterStringFontHeightF  = 0.018
res@tmXTLabelsOn  = False
res@tmYRLabelsOn  = False
vcres = res
res@tmXTOn        = False
res@tmYROn        = False
res@tmXBLabelFontHeightF = 0.023
res@tmYLLabelFontHeightF = 0.023
res@tmXBTickSpacingF = 45.0
res@tmYLTickSpacingF = 20.0




res@tmYLMajorThicknessF = 1.2
res@tmXBMajorThicknessF = 1.2
res@tmYLMinorThicknessF = 1.0
res@tmXBMinorThicknessF = 1.0
res@tmXBMajorLengthF = 0.018
res@tmYLMinorLengthF = 0.018
res@tmXBMinorLengthF = 0.016
res@tmYLMinorLengthF = 0.016
res@tmXBLabelDeltaF = -0.35
res@tmYLLabelDeltaF = -0.35
res@tmXBMajorOutwardLengthF = 0.018
res@tmYLMajorOutwardLengthF = 0.018
res@tmXBMinorOutwardLengthF = 0.016
res@tmYLMinorOutwardLengthF = 0.016





res@mpMinLatF = -90
res@mpMaxLatF = 90
res@mpMinLonF = 0
res@mpMaxLonF = 360
res@mpCenterLonF = 180
res@tmXBMode = "Explicit"
res@tmXBValues = (/-150,-120,-90,-60,-30,0,30,60,90,120,150,180/)
res@tmXBLabels = (/"150W","120W","90W","60W","30W","0","30E","60E","90E","120E","150E","180"/)
res@tmXBMinorOn = True
res@tmXBMinorValues = fspan(-180, 180+15, 25+1)

res@tmYLMinorOn = True
res@tmYLMinorValues = fspan(-90, 90, 19)

res@mpFillOn                    = False
res@mpOutlineOn                 = True
res@mpGeophysicalLineThicknessF = 1.5
res@mpGeophysicalLineColor      = "gray10";"gray20"
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
vcres@vcMinDistanceF          = 0.030             ; thin out vectors
vcres@vcLineArrowColor        = "black"           ; change vector color
vcres@vcRefAnnoOn             = True
vcres@vcRefLengthF            = 0.030            ; ref vec length
vcres@vcRefAnnoFontHeightF    = 0.019   ;参考箭头字体大小
vcres@vcRefAnnoString1On      = True
vcres@vcRefAnnoString2On      = False
vcres@vcLineArrowThicknessF   = 1.0            ; make vectors larger
vcres@vcVectorDrawOrder = "PostDraw"
vcres@vcRefAnnoOrthogonalPosF = -0.20 ;;正向下
vcres@vcRefAnnoParallelPosF   = 1.0  ;;正向右
vcres@vcRefAnnoString1        = "10"
vcres@vcRefMagnitudeF         = 10



sres = True
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
sres@cnLineColor = "red"
sres@cnLevelSelectionMode = "ExplicitLevels"
sres@cnLevels         = (/2700.0/)
sres@cnFillColors     = (/"red","gray20"/)
sres@cnMissingValFillColor= -1
sres@cnLineThicknessF = 3.0
sres@cnLineDashPattern = 14


sres2@cnLevelSelectionMode = "ExplicitLevels"
sres2@cnLevels         = (/2/)
sres2@cnLineColor = "blue"
sres2@cnLineDashPattern = 14
sres2@cnLineThicknessF = 4
; sres2@cnLineColors     = (/"red","red","red","red"/)
; sres2@cnLineDashPatterns = (/16,16,16,0,0,0,0/)
; sres2@cnLineThicknesses = (/4,4,4,4,4,4,4/)









;;*************************打点***************************;;
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
ores@cnMonoFillPattern    = False            ; want multiple patterns                                                               
; ores@cnFillPatterns       = (/-1,11/)     ; the patterns                                                                         
ores@cnMonoFillScale      = False            ; want different densities                                                             
; ores@cnFillScales         = (/0.6,0.6/)    ; change densities                                                                         
ores@cnMonoFillColor      =True                                                                                                     
ores@cnFillDotSizeF       = 0.003    
ores@cnFillDrawOrder ="postDraw"
ores@cnFillColor = "white"


ores@cnMonoFillPattern    = False            ; want multiple patterns                                                                                                                                     
ores@cnMonoFillScale      = False            ; want different densities   
ores@cnFillScales         = (/0.8,0.8,0.8/)    ; change densities  


ores1 = ores
ores1@cnFillColor = "black"
ores1@cnLevels        = (/0.05/)   ;;95%
ores1@cnFillPatterns  = (/3,-1/)     ; the patterns   
 
ores2 = ores
ores2@cnFillColor = "white"
ores2@cnLevels        = (/0.05/)   ;;95%, 0.294
ores2@cnFillPatterns  = (/3,-1/)     ; the patterns  

ores3 = ores
ores3@cnFillColor = "black"
ores3@cnLevels        = (/0.10/)   ;;95%, 0.294
ores3@cnFillPatterns  = (/-1,3/)     ; the patterns 

ores4 = ores
ores4@cnFillColor = "white"
ores4@cnLevels        = (/0.10/)   ;;95%, 0.294
ores4@cnFillPatterns  = (/-1,3/)     ; the patterns 


;**********************plot**************************
plot = new(14, "graphic")
vector = new(14, "graphic")
contour1 = new(14, "graphic")
contour2 = new(14, "graphic")
contour3 = new(14, "graphic")
contour4 = new(14, "graphic")
contour5 = new(14, "graphic")
contour6 = new(14, "graphic")
contour7 = new(14, "graphic")
contour8 = new(14, "graphic")
topooo = new(14, "graphic") 











;;---------------------------------------------------------------------------------------------------;;
pltType = "eps"
pltName = "/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure4_CMIP6_CESM"
wks = gsn_open_wks(pltType, pltName)


res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"
res@tmXBLabelsOn  = False
res@tmXBLabelsOn  = True
res@mpGeophysicalLineThicknessF = 0.8



res@mpMinLonF = 30
res@mpMaxLonF = 135
res@mpMinLatF = -10
res@mpMaxLatF = 50
res@mpCenterLonF = 0


res@cnFillMode = "CellFill"
res@cnFillPalette = "MPL_RdBu"
res@cnLevels     := (/-0.5,-0.3,-0.1,0.1,0.3,0.5/)
res@cnFillColors := (/100,88,77,-1,48,37,25/)
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels     := (/-0.5,-0.3,-0.1,0.1,0.3,0.5/)
res@cnFillColors := (/50,70,95,-1,160,185,205/)
plot(0) = gsn_csm_contour_map_ce(wks, TSMX_NIO_diff_avg, res)

delete(res@cnFillMode)
res@cnFillPalette = "MPL_RdBu"
res@cnLevels     := (/-1.0,-0.5,-0.1,0.1,0.5,1.0/)
res@cnFillColors := (/100,88,77,-1,48,37,25/)
res@cnFillColors := (/100,88,77,-1,46,37,29/)
plot(1) = gsn_csm_contour_map_ce(wks, wbt_NIO_diff_avg, res)
contour1(1) = gsn_csm_contour(wks, wbt_NIO_t, ores1)
overlay(plot(1), contour1(1))

res@cnFillPalette = "MPL_BrBG"
res@cnLevels     := (/-1.0,-0.5,-0.1,0.1,0.5,1.0/)
res@cnFillColors := (/28,37,50,-1,75,89,98/)
plot(2) = gsn_csm_contour_map_ce(wks, Q_NIO_diff_avg, res)
contour1(2) = gsn_csm_contour(wks, Q_NIO_t, ores2)
overlay(plot(2), contour1(2))
vcres@vcRefAnnoString1 = "12"
vcres@vcRefMagnitudeF  = 12
vector(2) = gsn_csm_vector(wks, QU_NIO_diff_avg, QV_NIO_diff_avg, vcres)
overlay(plot(2), vector(2))







; res@cnFillMode = "CellFill"
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevels     := (/-15,-10,-5,5,10,15/)
res@cnFillColors := (/50,70,95,-1,160,185,205/)
res@cnFillColors := (/55,75,95,-1,160,180,200/)
plot(3) = gsn_csm_contour_map_ce(wks, SA_regSST, res)
keyres = res
keyres@cnFillPalette = "GMT_gray"
keyres@cnLevels     := (/0.5/)
keyres@cnFillColors := (/-1,7/)
keyres@cnInfoLabelOn = False
contour1(3) = gsn_csm_contour(wks, key_region3, keyres)
overlay(plot(3), contour1(3))
contour2(3) = gsn_csm_contour(wks, SA_corSST, ores3)
overlay(plot(3), contour2(3))


delete(res@cnFillMode)
res@cnFillPalette = "MPL_RdBu"
res@cnLevels     := (/-30,-20,-10,10,20,30/)
res@cnFillColors := (/100,88,77,-1,48,37,27/)
plot(4) = gsn_csm_contour_map_ce(wks, SA_regwbt, res)
contour1(4) = gsn_csm_contour(wks, SA_corwbt, ores3)
overlay(plot(4), contour1(4))

res@cnFillPalette = "MPL_BrBG"
res@cnLevels     := (/-25,-15,-5,5,15,25/)
res@cnFillColors := (/28,37,50,-1,75,89,98/)
plot(5) = gsn_csm_contour_map_ce(wks, SA_regq, res)
contour1(5) = gsn_csm_contour(wks, SA_corq, ores4)
overlay(plot(5), contour1(5))
vcres@vcRefAnnoString1 = "3"
vcres@vcRefMagnitudeF  = 3
vcres@vcMinDistanceF   = 0.035
vector(5) = gsn_csm_vector(wks, SA_regqu, SA_regqv, vcres)
overlay(plot(5), vector(5))



txres               = True
txres@txFont = "helvetica-bold"
txres@txPerimOn = False
txres@txFontHeightF = 0.035
txres@txFontColor = "red4"
dtext1 = gsn_add_text(wks, plot(2), "C", 64, 18.5, txres)
dtext2 = gsn_add_text(wks, plot(5), "C", 62, 20.0, txres)





lnres   =    True
lnres@gsLineThicknessF = 2.5
lnres@gsLineDashPattern=0

; SAsia_lat1 = 40
; SAsia_lat2 = 65
; SAsia_lon1 = 120
; SAsia_lon2 = 140
; y1=SAsia_lat1
; y2=SAsia_lat2
; x1=SAsia_lon1
; x2=SAsia_lon2
; x=(/x1,x2,x2,x1,x1/)
; y=(/y2,y2,y1,y1,y2/)
; lnres@gsLineColor="green4"
;   dum11 = gsn_add_polyline(wks,plot(0),x,y,lnres)

; NA_lat1 = 30
; NA_lat2 = 40
; NA_lon1 = -105 + 360
; NA_lon2 = -85 + 360
; y1=NA_lat1
; y2=NA_lat2
; x1=NA_lon1
; x2=NA_lon2
; x=(/x1,x2,x2,x1,x1/)
; y=(/y2,y2,y1,y1,y2/)
; lnres@gsLineColor="green4"
;   dum12 = gsn_add_polyline(wks,plot(3),x,y,lnres)



pres = True
pres@gsnFrame = False
pres@gsnPanelLabelBar = True
pres@lbOrientation        = "Vertical"
pres@pmLabelBarWidthF     = 0.03
pres@pmLabelBarHeightF    = 0.13
; pres@pmLabelBarParallelPosF = 0.11
pres@pmLabelBarOrthogonalPosF = 0.004
pres@lbLabelFontHeightF  = 0.0087


pres@gsnPanelTop  = 0.715
pres@gsnPanelBottom = 0.465
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.34
pres@pmLabelBarOrthogonalPosF = 0.005
gsn_panel(wks,plot(3),(/1,1/),pres)
pres@gsnPanelLeft  = 0.34
pres@gsnPanelRight = 0.67
gsn_panel(wks,plot(4),(/1,1/),pres)
pres@gsnPanelLeft  = 0.67
pres@gsnPanelRight = 1.0
pres@pmLabelBarOrthogonalPosF = 0.007
gsn_panel(wks,plot(5),(/1,1/),pres)

pres@gsnPanelTop  = 0.445
pres@gsnPanelBottom = 0.195
pres@gsnPanelLeft  = 0.01
pres@gsnPanelRight = 0.34
pres@pmLabelBarOrthogonalPosF = 0.005
gsn_panel(wks,plot(0),(/1,1/),pres)
pres@gsnPanelLeft  = 0.34
pres@gsnPanelRight = 0.67
gsn_panel(wks,plot(1),(/1,1/),pres)
pres@gsnPanelLeft  = 0.67
pres@gsnPanelRight = 1.0
pres@pmLabelBarOrthogonalPosF = 0.007
gsn_panel(wks,plot(2),(/1,1/),pres)




txres               = True
txres@txFontHeightF = 0.013
txres@txFont = "helvetica"
txres@txFontColor = "black"
gsn_text_ndc(wks, "CESM    response    to    North    Indian    Ocean    warming",  0.52, 0.465, txres)
gsn_text_ndc(wks, "CMIP6    regression    onto    South    Asian    wet-bulb    temperature",  0.52, 0.740, txres)


txres@txFontHeightF = 0.0105
gsn_text_ndc(wks, "SST",  0.16, 0.685, txres)
gsn_text_ndc(wks, "Wet-bulb   temperature",  0.49, 0.685, txres)
gsn_text_ndc(wks, "q850   &   quv850",  0.82, 0.685, txres)
gsn_text_ndc(wks, "SST   forcing",  0.16, 0.415, txres)
gsn_text_ndc(wks, "Wet-bulb   temperature",  0.49, 0.415, txres)
gsn_text_ndc(wks, "q850   &   quv850",  0.82, 0.415, txres)



txres@txFont = "helvetica-bold"
txres@txFontColor = "black"
txres@txFontHeightF = 0.0145
txres@txAngleF = 0.0
gsn_text_ndc(wks, "a",  0.03, 0.70, txres)
gsn_text_ndc(wks, "b",  0.36, 0.70, txres)
gsn_text_ndc(wks, "c",  0.69, 0.70, txres)
gsn_text_ndc(wks, "d",  0.03, 0.43, txres)
gsn_text_ndc(wks, "e",  0.36, 0.43, txres)
gsn_text_ndc(wks, "f",  0.69, 0.43, txres)




frame(wks)
   delete(wks)  
   system("convert -geometry 2000x2000 -density 1600x1600 -trim " + pltName + "."+pltType + " " + pltName + ".png")
   system("rm " + pltName + "." + pltType + " -f")

end








