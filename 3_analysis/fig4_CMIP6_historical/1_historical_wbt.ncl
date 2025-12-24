load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin


f1 = addfile("/public/home/fcai/CMIP6/historical_tas_day/tas_for8humidheat/historical_tas_NH_JJA_30models_1979_2014.nc","r")
tas = doubletofloat(f1->tas)
tas@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(tas))) then
      replace_ieeenan (tas, value, 0)
      tas@_FillValue = value
end if
printVarSummary(tas)
printMinMax(tas, 1)  ;; 329 K

f2 = addfile("/public/home/fcai/CMIP6/historical_hurs_day/hurs_for8humidheat/historical_hurs_NH_JJA_30models_1979_2014.nc","r")
hurs = doubletofloat(f2->hurs)
hurs@_FillValue = 1.0e20
value = 1.0e20
if (any(isnan_ieee(hurs))) then
      replace_ieeenan (hurs, value, 0)
      hurs@_FillValue = value
end if
printVarSummary(hurs)
printMinMax(hurs, 1)



wbt = wetbulb_stull(tas, hurs, (/1,0/), False)
copy_VarCoords(tas, wbt)
printMinMax(wbt, 1)



;;-------------------- to nc file ---------------------;;
output="/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig4_CMIP6_historical/"
name1="wbt_historical_30models_1979_2014"
system("/bin/rm -f "+output+name1+".nc")
ncdf1=addfile(output+name1+".nc","c")
ncdf1->wbt=wbt


end

