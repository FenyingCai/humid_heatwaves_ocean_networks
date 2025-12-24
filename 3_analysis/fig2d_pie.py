

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt




##--##--##--##--##--##--       --##--##--##--##--##--##
ds1 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/Flow_landocean_lag0_5_NH.nc")

flow0 = ds1.flow_weight[0,:,:].data
NH_land1_day0 = np.zeros(4)
NH_land2_day0 = np.zeros(4)
for i in range(4): 
     NH_land1_day0[i] = flow0[2,i] / np.nansum(flow0[2,:], axis=0)
     NH_land2_day0[i] = flow0[3,i] / np.nansum(flow0[3,:], axis=0)

flow5 = ds1.flow_weight[1,:,:].data
NH_land1_day5 = np.zeros(4)
NH_land2_day5 = np.zeros(4)
for i in range(4): 
     NH_land1_day5[i] = flow5[2,i] / np.nansum(flow5[2,:], axis=0)
     NH_land2_day5[i] = flow5[3,i] / np.nansum(flow5[3,:], axis=0)





ds1 = xr.open_dataset("/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/fig2_flow/Flow_landocean_lag0_5_SH.nc")

flow0 = ds1.flow_weight[0,:,:].data
SH_land1_day0 = np.zeros(4)
SH_land2_day0 = np.zeros(4)
for i in range(4): 
     SH_land1_day0[i] = flow0[2,i] / np.nansum(flow0[2,:], axis=0)
     SH_land2_day0[i] = flow0[3,i] / np.nansum(flow0[3,:], axis=0)

flow5 = ds1.flow_weight[1,:,:].data
SH_land1_day5 = np.zeros(4)
SH_land2_day5 = np.zeros(4)
for i in range(4): 
     SH_land1_day5[i] = flow5[2,i] / np.nansum(flow5[2,:], axis=0)
     SH_land2_day5[i] = flow5[3,i] / np.nansum(flow5[3,:], axis=0)









# %%
## Plotting Import
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('agg')

## Define a Map Function for plotting
def create_map(ax):
  # plt.axes( xscale="log" )
  plt.xlim(-1.0,37.0)
  plt.ylim(2.0,5.5)
  ax.tick_params(axis='both', length=6.0, width=1.0, labelsize=8.5)
  ax.spines['bottom'].set_linewidth(1.0)
  ax.spines['top'].set_linewidth(1.0)
  ax.spines['left'].set_linewidth(1.0)
  ax.spines['right'].set_linewidth(1.0)
  ax.spines['left'].set_color('black')
  ax.tick_params(axis='y', right=False, length=6.0, width=1.0, labelsize=9.5, color='black')
  return ax


##--##--##--  create figure  --##--##--##
width_in_inches = 210 / 25.4  ## A4 page
height_in_inches = 297 / 25.4

fig = plt.figure(dpi=400,figsize=(width_in_inches,height_in_inches))

ax5 = plt.axes([0.03,0.30, 0.15,0.2])
ax5.pie(SH_land1_day0,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f5d33d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比'
        labeldistance=0.65, pctdistance=0.65, 
       )
ax6 = plt.axes([0.18,0.30, 0.15,0.2])
ax6.pie(SH_land1_day5,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f5d33d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比'
        labeldistance=0.65, pctdistance=0.65, 
       )

ax7 = plt.axes([0.36,0.30, 0.15,0.2])
ax7.pie(SH_land2_day0 ,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f5d33d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比
        labeldistance=0.7, pctdistance=0.65, 
       )
ax8 = plt.axes([0.51,0.30, 0.15,0.2])
ax8.pie(SH_land2_day5,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f5d33d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比
        labeldistance=0.7, pctdistance=0.65, 
       )





ax5.text(-1.1,1.9, 'f', fontsize=15.0, fontweight='bold', fontfamily='sans-serif')
ax5.text(0.05,1.7, 'SH   coastland', fontsize=12.0)
ax5.text(-0.30,1.3, 'humid   heatwaves', fontsize=12.0)
ax5.text(-0.82,-1.40, 'lag = 0 day', fontsize=10.5)
ax6.text(-0.80,-1.40, 'lag = 5 day', fontsize=10.5)
ax5.text(-0.35,-1.85, '54.3%', fontsize=10.5)
ax6.text(-0.35,-1.85, '67.2%', fontsize=10.5)


ax7.text(-1.1,1.9, 'g', fontsize=15.0, fontweight='bold', fontfamily='sans-serif')
ax7.text(0.30,1.7, 'SH   inland', fontsize=12.0)
ax7.text(-0.25,1.3, 'humid   heatwaves', fontsize=12.0)
ax7.text(-0.82,-1.40, 'lag = 0 day', fontsize=10.5)
ax8.text(-0.80,-1.40, 'lag = 5 day', fontsize=10.5)
ax7.text(-0.35,-1.85, '45.6%', fontsize=10.5)
ax8.text(-0.35,-1.85, '62.2%', fontsize=10.5)


fig.savefig('/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure2d_pie.png', bbox_inches = 'tight')
# fig.savefig('/public/home/fcai/extreme1_AT/NC2_1980_2010/3_era5_analysis/Figure1.pdf')
plt.show()
















