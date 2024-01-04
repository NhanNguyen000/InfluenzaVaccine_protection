# code from Saumya
load("/vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/HAIreclassify2.RData")
imed_metaData = HAIreclassify2$iMED
metaFile = read.csv("/vol/projects/CIIM/Influenza/iMED/proteomic/analysis/output/cohort2_metadata.csv",header = T)
metaFileMerged = merge(metaFile,imed_metaData,by.x="PatientID",by.y="patientID")

load("/vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")
metabols = mebo_Dat$iMED_2015
metaFileMerged = metaFileMerged[match(rownames(metabols),metaFileMerged$name),]
metabols = t(metabols)

metaFileMerged$H1N1_reclassify = factor(metaFileMerged$H1N1_reclassify,levels = c("HH","LH","HL","LL"))
metaFileMerged$H1N1_reclassifyTime = paste0(metaFileMerged$H1N1_reclassify,"_",metaFileMerged$time)

modTime = model.matrix(~0 + H1N1_reclassifyTime + age + gender.x, data = metaFileMerged)
corfit = duplicateCorrelation(metabols,modTime,block = metaFileMerged$PatientID)
corfit$consensus.correlation
fitTime = lmFit(metabols,modTime,block = metaFileMerged$PatientID,correlation = corfit$consensus.correlation)
fitTime = eBayes(fitTime)
summary(decideTests(fitTime))
cont.matrix = makeContrasts(HH_T2vsT1 = H1N1_reclassifyTimeHH_T2 - H1N1_reclassifyTimeHH_T1,
                            HH_T3vsT1 = H1N1_reclassifyTimeHH_T3 - H1N1_reclassifyTimeHH_T1,
                            LH_T2vsT1 = H1N1_reclassifyTimeLH_T2 - H1N1_reclassifyTimeLH_T1,
                            LH_T3vsT1 = H1N1_reclassifyTimeLH_T3 - H1N1_reclassifyTimeLH_T1,
                            HL_T2vsT1 = H1N1_reclassifyTimeHL_T2 - H1N1_reclassifyTimeHL_T1,
                            HL_T3vsT1 = H1N1_reclassifyTimeHL_T3 - H1N1_reclassifyTimeHL_T1,
                            LL_T2vsT1 = H1N1_reclassifyTimeLL_T2 - H1N1_reclassifyTimeLL_T1,
                            LL_T3vsT1 = H1N1_reclassifyTimeLL_T3 - H1N1_reclassifyTimeLL_T1, levels = modTime)

fit2 = contrasts.fit(fitTime,cont.matrix)
fit2 = eBayes(fit2)
summary(decideTests(fit2))

HH_T2vsT1_metabol = topTable(fit2,coef = 1,number = Inf)
HH_T3vsT1_metabol = topTable(fit2,coef = 2,number = Inf)
LH_T2vsT1_metabol = topTable(fit2,coef = 3,number = Inf)
LH_T3vsT1_metabol = topTable(fit2,coef = 4,number = Inf)
HL_T2vsT1_metabol = topTable(fit2,coef = 5,number = Inf)
HL_T3vsT1_metabol = topTable(fit2,coef = 6,number = Inf)
LL_T2vsT1_metabol = topTable(fit2,coef = 7,number = Inf)
LL_T3vsT1_metabol = topTable(fit2,coef = 8,number = Inf)

####### these are the numbers from summary(decideTests(fit2)) for plotting
down = c(0,-45,-44,0,-76,-100,0,-34,-55,0,-69,-70)
up =   c(0,13,23,0,39,66,0,37,36,0,34,50)
time = c("T1","T2","T3","T1","T2","T3","T1","T2","T3","T1","T2","T3")
Category = c("HH","HH","HH","LH","LH","LH","HL","HL","HL","LL","LL","LL")
newMat = data.frame(Down = down,Up = up, time = time,Category = Category)
newMatMelted = melt(newMat)
newMatMelted$CategoryVariable = paste0(newMatMelted$Category,"_",newMatMelted$variable)

ggplot(newMatMelted,aes(x=time,y=value,group=CategoryVariable))+
  geom_line(position="identity",aes(col=Category),size = 1)+
  geom_point(aes(col=Category),size=3)+theme_classic()+
  #ylim(-7,7)+#+scale_color_manual(values = c("dodgerblue","indianred2"))+
  geom_hline(yintercept = 0,linetype="dashed",color = "darkgray")+
  ylab("number of metabolites")