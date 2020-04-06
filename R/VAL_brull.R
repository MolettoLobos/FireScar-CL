library(rgdal)
library(raster)
library(ggplot2)
library(tidyverse)
library(caret)
rm(list=ls())


## load your polygone shapefile
#dsn<-"C:/Users/italo/Downloads/italo/italo/SHP/BD_incendios.shp"
#dsn = "C:/Users/italo/Downloads/italo/italo/SHP/fire_subset.shp"
#dsn = ('C:/Users/italo/Documents/MolettoLobos/LEPFOR_documentos/PerimetrosIncendiosCONAF26092018/2015-2016/Perimetros_IFM_CONAF_2015-2017_Attr.shp')
dsn = ('C:/Users/italo/Google Drive/FireScars_KML_indrnbr5/output/Perimetros_IFM_CONAF_2015-2017_Attr.shp')
#ogrListLayers(dsn) #list the layers in the above directory
shp_truth_raw <-readOGR(dsn) #shapefile name is lol
shp_truth_raw$layer = shp_truth_raw$FireScar16
df_truth_raw = shp_truth_raw@data[!duplicated(shp_truth_raw$layer),]

master_filter_raw=read.csv('C:/Users/italo/Google Drive/FireScars_KML_indrnbr5/output/data_incendios.csv',sep=";")
master_filter <- filter(master_filter_raw, 
                        master_filter_raw$cumple %in% c("si","no y si"))

df_truth <- filter(df_truth_raw, 
                   df_truth_raw$Incendio %in% master_filter$id)
df_truth = df_truth[!duplicated(df_truth$Incendio),]

df_evaluation = merge(master_filter,df_truth,by.x="id",by.y="Incendio",all=TRUE)
df_evaluation<-df_evaluation[!(df_evaluation$id=="109 - TERMA INTERNACIONAL"),]
df_truth_raw_dates=merge(x = shp_truth_raw, y = master_filter, by.x = "Incendio",by.y="id", all.x = TRUE)
df_truth_filter_nadates = df_truth_raw_dates[!is.na(df_truth_raw_dates$fecha_inicio), ]
#check your data
#original path
#path = "C:/Users/italo/Google Drive/FireScars_KML_val_original"
#path = "C:/Users/italo/Google Drive/FireScars_KML_indrnbr5"
path = 'C:/Users/italo/Google Drive'
#path = 'C:/Users/italo/Documents/'
path = 'C:/Users/italo/Downloads/Google engine-20200129T122211Z-001'

cal_list = list.files(path = path, patter = "FireScars_KML", full.names = T)
cal_list = list.files(path = path, patter = 'Google', full.names = T)
#cal_list = list.files(path = path, patter = "FireScars_KML", full.names = T)
n = length(cal_list)
df_matrix = data.frame(r2 = vector(mode="double", length=n), RMSE=vector(mode="double", length=n),BIAS = vector(mode="double",length = n),SIGMA = vector(mode="double",length = n),n = vector(mode="double",length = n))
for (run_kml in 1:n){
cal_name = strsplit(cal_list[run_kml],'/')[[1]][5]
#validation path
#threshold = 'threshold_nbr_011'
#path = paste0("C:/Users/italo/Google Drive/FireScars_KML_val_",threshold)
shp_model_list = list.files(path = cal_list[run_kml], patter = ".kml", full.names = T)
x=length(df_evaluation$id)
df = data.frame(area = vector(mode="double", length=x), ID=vector(mode="numeric", length=x),area_truth = vector(mode="double",length = x), ID_Brull=vector(mode="numeric", length=x),no_fire_pix=vector(mode="numeric", length=x),fire_pix=vector(mode="numeric", length=x),comission=vector(mode="numeric", length=x),omission=vector(mode="numeric", length=x),accuracy=vector(mode="numeric", length=x),kappa=vector(mode="numeric", length=x))
shp_truth_raw$id = as.character(shp_truth_raw$Incendio)
#remover data no coincidente extra

for (i in 1:x) {
  #id 77 = LEVICAN
  id_site=as.character(df_evaluation$id[i])
  date_site=as.character(df_evaluation$fecha_inicio[i])
  
  shp_truth_raw_idsubset=subset(df_truth_filter_nadates,df_truth_filter_nadates$Incendio==id_site)
  shp_truth_raw_idsubset_datesubset=subset(shp_truth_raw_idsubset,shp_truth_raw_idsubset$Incendio==id_site)
  #firescars_area = shp_truth_raw_idsubset_datesubset[!duplicated(shp_truth_raw_idsubset_datesubset$layer),]
  firescars_area = shp_truth_raw[shp_truth_raw@data$Incendio == id_site,]
  poly_filter = paste0(as.character(firescars_area$layer),'.kml')
  #date_filter=paste0(strsplit(as.character(df_evaluation$layer[i]),'_')[[1]][5],'.kml')
  kml_file=list.files(path = cal_list[run_kml], patter = poly_filter, full.names = T)[1]
  if (poly_filter == '.kml') {next}
  #print(length(dummy))}
  #substring(shp_model_list,77,100)
  #model
  #df_model[i,2] = as.numeric(substring(strsplit(shp_model_list[i],'_')[[1]][4],3,1000))
  #df_model[i,2] = as.numeric(substring(strsplit(shp_model_list[i],'_')[[1]][6],3,1000))
  #df_model[i,2] = as.numeric(substring(strsplit(shp_model_list[i],'_')[[1]][8],3,1000)) # points conaf
  df
  #df[i,2] = (substring(strsplit(kml_file,'_')[[1]][8],1,1000)) #brull
  df[i,2] = (substring(strsplit(kml_file,'_')[[1]][3],1,1000))
  #date_kml=as.numeric(substring(strsplit(kml_file,'_')[[1]][10],1,4))
  date_kml = as.numeric(substring(strsplit(kml_file,'_')[[1]][5],1,4))
  if (is.na(date_kml)) {
    print('date kml not found')
    next
  }
  if (date_kml < 2015) {
    print(date_kml)
    print('date kml not valid for evaluation')
    next
  }
  try(
    df[i,1] <-sum(area(readOGR(kml_file)))/10000
  ,TRUE)
  df[i,3] = as.numeric(as.character(firescars_area$sup_GIS))[1]
  id_inc = as.character(firescars_area$layer)[1]
  if (id_inc == 'Fire_scars_ID189238_Date_2017-2-21') {
    next #ZEMITA NO CUMPLE CON CONTROL DE CALIDAD
  }
  df[i,4] = id_site
  print(i)
  
  
  #CONFUSION MATRIX EVALUATION
  firescars_area$NUMBER=1
  sp_model=readOGR(kml_file)
  sp_model$NUMBER = 1
  projection_fire = CRS("+init=epsg:32719")
  #give metric to the file
  shp_truth_UTM19S <- spTransform(firescars_area,
                                  crs(projection_fire))
  shp_model_UTM19S <- spTransform(sp_model,
                                 crs(projection_fire))
  extent_fire_truth = extent(shp_truth_UTM19S)
  extent_fire_model = extent(shp_model_UTM19S)
  extent_fire = extent( min(c(extent_fire_truth@xmin,extent_fire_model@xmin)),max(c(extent_fire_truth@xmax,extent_fire_model@xmax)),min(c(extent_fire_truth@ymin,extent_fire_model@ymin)), max(c(extent_fire_truth@ymax,extent_fire_model@ymax))     )
  print(extent_fire)
  resolution_fire = 30 #meters
  nrows_fire = as.integer((extent_fire@ymax-extent_fire@ymin)/resolution_fire)
  ncols_fire = as.integer((extent_fire@xmax-extent_fire@xmin)/resolution_fire)
  ## raster template
  rst_template <- raster(ncols = ncols_fire, nrows = nrows_fire, 
                         crs = projection(projection_fire), 
                         ext = extent(extent_fire))
  rst_truth <- rasterize(shp_truth_UTM19S, rst_template,field=shp_truth_UTM19S$NUMBER)
  rst_truth[is.na(rst_truth[])] <- 0 
  rst_model = rasterize(shp_model_UTM19S, rst_template,field=shp_model_UTM19S$NUMBER)
  rst_model[is.na(rst_model[])] <- 0 
  rst_evaltruth = rst_truth+10
  crosstab_fire = rst_model+rst_evaltruth 
  crosstab_fire[rst_model == 0 & rst_evaltruth == 11 ] = 13 #error de omision
  crosstab_fire_shp = rasterToPolygons(crosstab_fire,dissolve=TRUE)
  crosstab_fire_shp$ID_Brull = id_site
  crosstab_fire_shp$ID = id_inc
  crosstab_fire_shp$Val = ""
  crosstab_fire_shp$Val[crosstab_fire_shp$layer == 10] = 'No incendio'
  crosstab_fire_shp$Val[crosstab_fire_shp$layer == 11] = 'Error de comision'
  crosstab_fire_shp$Val[crosstab_fire_shp$layer == 12] = 'Incendio'
  crosstab_fire_shp$Val[crosstab_fire_shp$layer == 13] = 'Error de omision'
  dsn_export = paste0("C:/Users/italo/Documents/MolettoLobos/LEPFOR_documentos/Resultados_incendio/Shp/crosstab_shape_",df[i,2],'.shp')
  writeOGR(obj=crosstab_fire_shp, dsn=dsn_export, layer="torn", driver="ESRI Shapefile") # this is in geographical projection
  
  #get confusion matrix
  model_df=as.data.frame(rst_model)
  truth_df = as.data.frame(rst_truth)
  model=factor(model_df$layer)
  truth = factor(truth_df$layer)
  eval = confusionMatrix(model,truth )
  conf_mat = eval[2]$table
  df[i,6]=conf_mat[4] #Pixel de incendio bien identificado
  df[i,5]=conf_mat[1] #Pixel no incendio bien identificado
  df[i,8]=conf_mat[3] #Error de Omisión
  df[i,7]=conf_mat[2] #Error de Comision
  df[i,9]=eval$overall[1]
  df[i,10]=eval$overall[2]
  }

#df_truth = data.frame(are  a=shp_truth$sup_GIS,KEY=shp_truth$Incendio)
#df_truth$KEY = gsub('\\.','',as.character(df_truth$KEY))
#df_model$ID = gsub('\\Ñ','N',as.character(df_model$ID))
#df = merge(df_truth,df_model, by.x="KEY", by.y="ID")
#df$area.x[df$area.x == 0] = NA
#LIMPIEZA DE DATOS
#eliminar key ID 177196 ya que se considera como punto NO VALIDO

df$res = df$area-df$area_truth
##Go through each row and determine if a value is zero
row_sub = apply(df, 1, function(row) all(row !=0 ))
##Subset as usual
df=df[row_sub,]


bias = mean(df$res,na.rm=TRUE)
sigma = sd(df$res,na.rm=TRUE)
rmse = sqrt(bias^2+sigma^2)

strRMSE = as.character(round(rmse,2))
strBIAS = as.character(round(bias,2))
strSIGMA = as.character(round(sigma,2))
lmodel = lm(area_truth~area,df)
int = format(coefficients(lmodel)[1],digits = 2)
slope = format(coefficients(lmodel[2]),digits = 2)
r2 = format(summary(lmodel)$r.squared, digits = 2)
df_val = data.frame(area_y=df$area,area_x=df$area_truth)
#View(df)

df_matrix[run_kml,1] = r2
df_matrix[run_kml,2] = rmse
df_matrix[run_kml,3] = bias
df_matrix[run_kml,4] = sigma
df_matrix[run_kml,5] = x

df$KEY[df$res == max(df$res,na.rm=T)]

lm_eqn <- function(df_val){
  m <- lm(area_y ~ area_x, df_val);
  eq <- substitute(italic(y) == a +b*italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = as.character(format(coef(m)[1], digits = 2)), 
                        b = as.character(format(coef(m)[2], digits = 2)), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
  
  graph=ggplot(data = df_val, aes(x = area_x, y = area_y)) + 
    geom_point(color='red',size=3) +
    geom_smooth(method = "lm", se = FALSE) +
    ggtitle('Incendio CONAF VS Landsat 8') +
    xlab('Area incendio CONAF [ha]') + ylab('Area incendio Landsat [ha]')+
    xlim(0,7500) + ylim(0,7500) +
    #expand_limits(x = 0, y = 0) +
    #scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    geom_text (x=5000,y=6000,label=lm_eqn(df_val),parse=TRUE,size=12)+
    annotate("text",x=5000,y=5500,label=paste0("RMSE = ",strRMSE),size=12) +
    annotate("text",x=5000,y=5000,label=paste0("SIGMA = ",strSIGMA),size=12) +
    annotate("text",x=5000,y=4500,label=paste0("BIAS = ",strBIAS),size=12) +
    theme_bw() +theme(axis.text=element_text(size=20),plot.title = element_text(size = 40, face = "bold"),
                        axis.title=element_text(size=20,face="bold"))
  graph
  ggsave(plot=graph,width= 10,height= 10, filename=paste0('C:/Users/italo/Documents/MolettoLobos/LEPFOR_documentos/PPT/Graph_',cal_name,'.png' ))
  write.csv(df,paste0('C:/Users/italo/Documents/MolettoLobos/LEPFOR_documentos/PPT/df_ddbb',cal_name,'.csv'))
  print(run_kml)
}

write.csv(df_matrix,'C:/Users/italo/Documents/MolettoLobos/LEPFOR_documentos/PPT/df_eval.csv')
write.csv(df_truth_filter_nadates,'C:/Users/italo/Documents/MolettoLobos/LEPFOR_documentos/PPT/df_rawdata.csv')
