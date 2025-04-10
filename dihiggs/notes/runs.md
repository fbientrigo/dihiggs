# Run and Time data
Allows to know the $\mathcal O(N)$ in order to predict efficiency of parallel version

ftrigo@nasanote:~/dihiggs/app$ ./ParamScanGen_pll conf/coarse_parametros3.conf out/out_5pll.csv
Total iterations: 105000
15.22 min 

ftrigo@nasanote:~/dihiggs/app$ ./ParamScanGen_pll conf/coarse_parametros2.conf out/out_4pll.csv
Total iterations: 56700
7.43 min 

ftrigo@nasanote:~/dihiggs/app$ ./ParamScanGen_pll conf/coarse_parametros.conf out/out_3pll.csv
Total iterations: 126360
17.39 min 

ftrigo@nasanote:~/dihiggs/app$ ./ParamScanGen_pll conf/0129a.conf out/0129a.csv
Total iterations: 90720
11.46 min

ftrigo@nasanote:~/dihiggs/app$ ./ParamScanGen_pll conf/0129b.conf out/0129b.csv
Total iterations: 921600
Elapsed: 142.01 min


# Para el caso de Phys

Nota: parece haber un bug en como se calcula el tiempo total, dado el step_0

722400
Progress: 0.66% | Elapsed: 58.09 min 

237600
Progress: 1.52% | Elapsed: 16.53 min

