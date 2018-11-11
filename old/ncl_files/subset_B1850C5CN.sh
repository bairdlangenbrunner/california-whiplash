file_list=(
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.04020101-04991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.05000101-05991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.06000101-06991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.07000101-07991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.08000101-08991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.09000101-09991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.10000101-10991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.11000101-11991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.12000101-12991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.13000101-13991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.14000101-14991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.15000101-15991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.16000101-16991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.17000101-17991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.18000101-18991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.19000101-19991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.20000101-20991231.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.21000101-22001231.nc
)
saveas_list=(
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.04020101-04991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.05000101-05991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.06000101-06991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.07000101-07991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.08000101-08991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.09000101-09991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.10000101-10991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.11000101-11991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.12000101-12991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.13000101-13991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.14000101-14991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.15000101-15991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.16000101-16991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.17000101-17991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.18000101-18991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.19000101-19991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.20000101-20991231_CA_REGION.nc
b.e11.B1850C5CN.f09_g16.005.cam.h1.PRECT.21000101-22001231_CA_REGION.nc
)

#for i in ${hist_list[*]}; do
for ((i=0;i<${#file_list[@]};++i)); do
	ncea -d lat,25.,50. -d lon,220.,250. ${file_list[i]} ${saveas_list[i]}
	echo $i
done