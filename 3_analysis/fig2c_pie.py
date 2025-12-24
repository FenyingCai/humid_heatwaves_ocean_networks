

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

ax1 = plt.axes([0.03,0.5, 0.15,0.2])
ax1.pie(NH_land1_day0,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f3c80d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比'
        labeldistance=0.65, pctdistance=0.65, 
       )
ax2 = plt.axes([0.18,0.5, 0.15,0.2])
ax2.pie(NH_land1_day5,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f3c80d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比'
        labeldistance=0.65, pctdistance=0.65, 
       )

ax3 = plt.axes([0.36,0.5, 0.15,0.2])
ax3.pie(NH_land2_day0 ,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f3c80d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比
        labeldistance=0.7, pctdistance=0.65, 
       )
ax4 = plt.axes([0.51,0.5, 0.15,0.2])
ax4.pie(NH_land2_day5,
        labels=[' ',' ',' ',' '], # 设置饼图标签
        colors=["#2b2a76b3", "#73aa43b3", "#f3c80d", "#a21d2fB3"], # 设置饼图颜色
        autopct='%.1f%%', # 格式化输出百分比
        labeldistance=0.7, pctdistance=0.65, 
       )




ax1.text(-1.1,1.9, 'd', fontsize=15.0, fontweight='bold', fontfamily='sans-serif')
ax1.text(0.05,1.7, 'NH   coastland', fontsize=12.0)
ax1.text(-0.30,1.3, 'humid   heatwaves', fontsize=12.0)
ax1.text(-0.82,-1.40, 'lag = 0 day', fontsize=10.5)
ax2.text(-0.80,-1.40, 'lag = 5 day', fontsize=10.5)
ax1.text(-0.35,-1.85, '59.1%', fontsize=10.5)
ax2.text(-0.35,-1.85, '68.5%', fontsize=10.5)

ax3.text(-1.1,1.9, 'e', fontsize=15.0, fontweight='bold', fontfamily='sans-serif')
ax3.text(0.30,1.7, 'NH   inland', fontsize=12.0)
ax3.text(-0.25,1.3, 'humid   heatwaves', fontsize=12.0)
ax3.text(-0.82,-1.40, 'lag = 0 day', fontsize=10.5)
ax4.text(-0.80,-1.40, 'lag = 5 day', fontsize=10.5)
ax3.text(-0.35,-1.85, '49.7%', fontsize=10.5)
ax4.text(-0.35,-1.85, '63.6%', fontsize=10.5)


ax3.text(-0.65,2.80, 'Attribution   of   land-area   humid   heatwaves', fontsize=14.5)




fig.savefig('/public/home/fcai/extreme5_degree/25_humid_heat/3_analysis1982_2023/Figure2c_pie.png', bbox_inches = 'tight')
# fig.savefig('/public/home/fcai/extreme1_AT/NC2_1980_2010/3_era5_analysis/Figure1.pdf')
plt.show()


