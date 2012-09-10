sfsegyread tape=data/sigsbee/sigsbee2a_migvel.sgy tfile=migvel-t.rsf hfile=migvel-h bfile=migvel-b > migvel-raw.rsf
< migvel-raw.rsf sfscale rscale=0.001 | sfscale rscale=0.3048 | sfput o1=0 d1=0.00762 label1=z unit1=km o2=3.05562 d2=0.01143 label2=x unit2=km > migvel.rsf
sfsegyread tape=data/sigsbee/sigsbee2a_stratigraphy.sgy tfile=strvel-t.rsf hfile=strvel-h bfile=strvel-b > strvel-raw.rsf
< strvel-raw.rsf sfscale rscale=0.001 | sfscale rscale=0.3048 | sfput o1=0 d1=0.00762 label1=z unit1=km o2=3.048 d2=0.00762 label2=x unit2=km > strvel.rsf
sfsegyread tape=data/sigsbee/sigsbee2a_nfs.sgy tfile=data-t.rsf hfile=data-h bfile=data-b > data.rsf
< data-t.rsf sfdd type=float | sfheadermath output="10925+fldr*150" | sfwindow > shot_data-ss.rsf
< data-t.rsf sfdd type=float | sfheadermath output="offset" | sfwindow > shot_data-oo.rsf
< shot_data-ss.rsf sfmath output=input/150 > shot_data-si.rsf
< shot_data-oo.rsf sfmath output=input/75 > shot_data-oi.rsf
< shot_data-oi.rsf sfcat axis=2 space=n shot_data-si.rsf | sftransp | sfdd type=int > shot_data-os.rsf
< data.rsf sfintbin head=shot_data-os.rsf xkey=0 ykey=1 | sfput d2=0.075 d3=0.150 o3=10.95 label1=t label2=o label3=s > shot_data-raw.rsf
< shot_data-raw.rsf sfmutter half=false t0=1.0 v0=6 | sfput label1=t unit1=s o2=0 d2=0.02286 unit2=km o3=3.33756 d3=0.04572 unit3=km > shot_data.rsf
< migvel.rsf sfmask min=4.499 | sfdd type=float > smask.rsf
< smask.rsf sfgrey parallel2=n labelrot=n wantaxis=y title="" pclip=100 min1=0 max1=9.144 label1=z unit1=km min2=3.048 max2=27.432 label2=x unit2=km screenratio=0.375 screenht=5.10938 wantscalebar=n parallel2=n labelsz=6 labelfat=3 titlesz=12 titlefat=3 allpos=y > Fig/smask.vpl
< migvel.rsf sfmask max=1.5 | sfdd type=float > wmask.rsf
< wmask.rsf sfgrey parallel2=n labelrot=n wantaxis=y title="" pclip=100 min1=0 max1=9.144 label1=z unit1=km min2=3.048 max2=27.432 label2=x unit2=km screenratio=0.375 screenht=5.10938 wantscalebar=n parallel2=n labelsz=6 labelfat=3 titlesz=12 titlefat=3 allpos=y > Fig/wmask.vpl
< smask.rsf sfadd wmask.rsf | sfmath output="1-input" > lmask.rsf
< lmask.rsf sfgrey parallel2=n labelrot=n wantaxis=y title="" pclip=100 min1=0 max1=9.144 label1=z unit1=km min2=3.048 max2=27.432 label2=x unit2=km screenratio=0.375 screenht=5.10938 wantscalebar=n parallel2=n labelsz=6 labelfat=3 titlesz=12 titlefat=3 allpos=y > Fig/lmask.vpl
< smask.rsf sfmath output="exp(-((x1-3)^2 + (x2-10)^2)/2)" | sfmask min=0.5 | sfdd type=float > BOS.rsf
< migvel.rsf sfmath output="1" > den.rsf
sfspike nsp=1 mag=1 n1=1500 d1=0.008 o1=0 k1=13 | sfpad end1=1500 | sfricker1 frequency=10 | sfwindow n1=1500 | sfscale axis=123 | sfput label1=t > wav_tmp.rsf
< wav_tmp.rsf sfsinc d1=0.001 n1=9000 o1=0 | sfput o1=-0.1 > wavx.rsf
< wavx.rsf sftransp plane=12 > wav.rsf
sfspike nsp=1 mag=1 n2=348 d2=0.02286 o2=0 k2=40 l2=308 n1=9000 d1=0.001 o1=0.0 | sfsmooth rect2=20 repeat=4 | sftransp > dmask.rsf
< migvel.rsf sfwindow n2=1 f2=0 | sfwindow f1=379 | sfvel1d vel=1.5 mask=wmask.rsf > vel1d.rsf
< migvel.rsf sfwindow n2=1 f2=0 | sfwindow f1=379 | sfput o1=0.0 | sfmath output="input -x1*0.05" | sfvel1d vel=1.5 mask=wmask.rsf > vel1d_80.rsf
< shot_data.rsf sfwindow n3=1 f3=123 | sfsinc pattern=wavx.rsf | sfbandpass fhi=20 | sftransp | sfreverse axis=1 opt=i | sfadd mode=p dmask.rsf > data-123.rsf
< data-123.rsf sftransp plane=13 | sfrtoc | sfmath output="imag(0.0)+real(x1)" | sfdd type=float | sfwindow n2=1 n3=1 | sfwindow squeeze=y > ss-123.rsf
< data-123.rsf sftransp plane=23 | sfrtoc | sfmath output="imag(0.0)+real(x1+x2)" | sfdd type=float | sfwindow n3=1 | sfput n1=2 n2=348 > rr-123.rsf
< wav.rsf sfawefd2d ompchunk=1 ompnth=0 verb=y free=n snap=y jsnap=100 nb=250 dabc=y vel=vel1d_80.rsf den=den.rsf sou=ss-123.rsf rec=rr-123.rsf wfl=wts-123.rsf jsnap=10 oqx=8 > dd-123.rsf
< data-123.rsf sfawefd2d ompchunk=1 ompnth=0 verb=y free=n snap=y jsnap=100 nb=250 dabc=y vel=vel1d_80.rsf den=den.rsf sou=rr-123.rsf rec=rr-123.rsf wfl=wtr_rev-123.rsf jsnap=10 oqx=8 > dr-123.rsf
< wtr_rev-123.rsf sfreverse axis=3 opt=i > wtr-123.rsf
< wts-123.rsf sfxcor2d uu=wtr-123.rsf axis=3 verb=y ompnth=0 nbuf=100 > Img-123.rsf

