

ald = load('/home/jason/ald-lj/data_11x.dat');

str.nmd = '/home/jason/disorder/lj/alloy/10K/0.0/10x/NMD/1/work';
nmd(1).nmd=load(strcat(str.nmd,'/NMDdata.mat'));
nmd(1).sed=load(strcat(str.nmd,'/SEDdata.mat'));
nmd(1).sed = nmd_convert_data(nmd(1).nmd,nmd(1).sed);

h = figure;
loglog(...
    ald(:,1)*nmd(1).nmd.LJ.tau/nmd(1).nmd.constant.s2ps,...
    ald(:,2),'.',...
    nmd(1).sed.freq,nmd(1).sed.life,'.'...
    )
print(h,'-dpng','m_ald_lj_nmd_compare.png')