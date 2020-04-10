# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 15:13:27 2019

@author: italo
"""
import ee
import datetime
ee.Initialize()

#--------------------------------------------------------------------------------------------------------------------------------
#--Este script carga la base de datos de incendios de donde se obtiene la fecha de inicio y termino del incendio, ---------------
#--su ubicación estimada y su superficie estimada. Selecciona las imagenes pre y post incendio mediante una reducción de media,--
#--Crea un área de analisis en donde compara el pre y post indice NBR.-----------------------------------------------------------
#--Genera el perimetro de cada incendio y su severidad.-------------------------------------------------------------------------- 
#--------------------------------------------------------------------------------------------------------------------------------

###############/--- BOX 0 ---###########################/ 
#################################################/
#################################################/
## Definir funciones para filtrar las imagenes a utiliziar - Seleccion de datos satelitales ###
#################################################/
#################################################/

# Generar mascara de nubes y sombra de nubes
# Obtención del valor de los bits correspondiente a los valores por pixel de la calidad de las imagenes
def getQABits(image, start, end, newName):
    # Calculo de bits a extraer
   pattern = 0
   for i in range(start, end + 1):
        pattern += 2 ** i
        #Generación de banda de calidad por pixel en bit
        #nuevo nombre a imagen
        return image.select([0], [newName]).bitwiseAnd(pattern).rightShift(start)

# Función para enmascarar sombra de nubes y nieve (7,8 sombra de nube y 9,10 nieve y hielo)
def cloud_shadows(image): 
    #Seleccionar la banda de calidad
    QA = image.select(['BQA']);
    #Aplicar función 
    return getQABits(QA, 7,8, 'Cloud_shadows').eq(1)

# Función para enmascarar nubes
def clouds(image):
    #Seleccionar la banda de calidad
    QA = image.select(['BQA']);
    #Aplicar función 
    return getQABits(QA, 4,4, 'Cloud').eq(0);
                            
# Función para enmascarar nieve
def snow(image):
    #Seleccionar la banda de calidad
    QA = image.select(['BQA']);
    #Aplicar función 
    return getQABits(QA, 9,10, 'Snow').eq(1)
                          
#Generación función de mascara
def maskClouds(image):
    cs = cloud_shadows(image);
    c = clouds(image);
    s = snow(image);
    image = image.updateMask(cs);
    image = image.updateMask(s);
    return image.updateMask(c);
                      
# Agregar banda representando la fecha de adquisicion #https:#developers.google.com/earth-engine/ic_composite_mosaic
def adddate(image):
    return image.addBands(image.metadata('system:time_start'));


def get_INDEX_L8(img):
    image = img
    ndvi = ee.Image(image.expression('(NIR - Red) / (NIR + Red)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4')})).rename('NDVI')
    nbr = ee.Image(image.expression('(NIR - SWIR2) / (NIR + SWIR2)', {
        'NIR': image.select('B5'),
        'SWIR2': image.select('B7')})).rename('NBR')
    nbr2 = ee.Image(image.expression('(SWIR1 - SWIR2) / (SWIR1 + SWIR2)', {
        'SWIR1': image.select('B6'),
        'SWIR2': image.select('B7')})).rename('NBR2')
    bai = ee.Image(image.expression('1/((0.1-Red)**2+(0.06-NIR)**2)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4')})).rename('BAI')
    gemi = ee.Image(image.expression('((2*((NIR**2)-(Red**2)))+(1.5*NIR)+(0.5*Red))/(NIR+Red+0.5)*(1-0.25*((2*((NIR**2)-(Red**2)))+(1.5*NIR)+(0.5*Red))/(NIR+Red+0.5))-((Red-0.125)/(1-Red))', {
        'NIR': image.select('B5'),
        'Red': image.select('B4')})).rename('GEMI')
    mirbi = ee.Image(image.expression('10*SWIR2-9.8*SWIR1+2', {
        'SWIR1': image.select('B6'),
        'SWIR2': image.select('B7')})).rename('MIRBI')
    ndmi = ee.Image(image.expression('(NIR - SWIR1) / (NIR + SWIR1)', {
        'NIR': image.select('B5'),
        'SWIR1': image.select('B6')})).rename('NDMI')
    savi = ee.Image(image.expression('((NIR - Red) / (NIR + Red + L)) * (1 + L)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4'),
        'L': 0.5 })).rename('SAVI')
    evi = ee.Image(image.expression('2.5 * (NIR - Red) / ((NIR + C1 * Red - C2 * Blue) + L)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4'),
        'C1': 6,
        'C2': 7.5,
        'Blue': image.select('B2'),
        'L': 1}
                                    )).rename('EVI')

    gndvi = ee.Image(image.expression('(NIR - green) / (NIR + green)', {
        'NIR': image.select('B5'),
        'green': image.select('B3')}
                                      )).rename('GNDVI')

    msavi2 = ee.Image(image.expression('(2 * NIR + 1 - ((2 * NIR + 1)**2 - 8 * (NIR -Red))**(1/2))/2', {
        'NIR': image.select('B5'),
        'Red': image.select('B4')
    })).rename('MSAVI2')

    albedo = image.expression('(0.356*B+0.130*R+0.373*NIR+0.085*SW1+0.072*SW2-0.018)',
                              {'B': image.select('B2'), 'R': image.select('B4'),
                               'NIR': image.select('B5'), 'SW1': image.select('B6'),
                               'SW2': image.select('B7')}).rename('alfa')

    return image.addBands(ndvi).addBands(nbr).addBands(nbr2).addBands(bai).addBands(gemi).addBands(mirbi).addBands(ndmi).addBands(savi).addBands(evi).addBands(gndvi).addBands(msavi2).addBands(albedo).copyProperties(img, ['system:time_start', 'system:time_end','LANDSAT_ID'])


def get_INDEX_L57(img):
        image = img
        ndvi = ee.Image(image.expression('(NIR - Red) / (NIR + Red)', {
            'NIR': image.select('B4'),
            'Red': image.select('B3')})).rename('NDVI')
        nbr = ee.Image(image.expression('(NIR - SWIR2) / (NIR + SWIR2)', {
            'NIR': image.select('B4'),
            'SWIR2': image.select('B7')})).rename('NBR')
        nbr2 = ee.Image(image.expression('(SWIR1 - SWIR2) / (SWIR1 + SWIR2)', {
            'SWIR1': image.select('B5'),
            'SWIR2': image.select('B7')})).rename('NBR2')
        bai = ee.Image(image.expression('1/((0.1-Red)**2+(0.06-NIR)**2)', {
            'NIR': image.select('B4'),
            'Red': image.select('B3')})).rename('BAI')
        gemi = ee.Image(image.expression(
            '((2*((NIR**2)-(Red**2)))+(1.5*NIR)+(0.5*Red))/(NIR+Red+0.5)*(1-0.25*((2*((NIR**2)-(Red**2)))+(1.5*NIR)+(0.5*Red))/(NIR+Red+0.5))-((Red-0.125)/(1-Red))',
            {
                'NIR': image.select('B4'),
                'Red': image.select('B3')})).rename('GEMI')
        mirbi = ee.Image(image.expression('10*SWIR2-9.8*SWIR1+2', {
            'SWIR1': image.select('B5'),
            'SWIR2': image.select('B7')})).rename('MIRBI')
        ndmi = ee.Image(image.expression('(NIR - SWIR1) / (NIR + SWIR1)', {
            'NIR': image.select('B4'),
            'SWIR1': image.select('B5')})).rename('NDMI')
        savi = ee.Image(image.expression('((NIR - Red) / (NIR + Red + L)) * (1 + L)', {
            'NIR': image.select('B4'),
            'Red': image.select('B3'),
            'L': 0.5})).rename('SAVI')
        evi = ee.Image(image.expression('2.5 * (NIR - Red) / ((NIR + C1 * Red - C2 * Blue) + L)', {
            'NIR': image.select('B4'),
            'Red': image.select('B3'),
            'C1': 6,
            'C2': 7.5,
            'Blue': image.select('B1'),
            'L': 1}
                                        )).rename('EVI')
        gndvi = ee.Image(image.expression('(NIR - green) / (NIR + green)', {
            'NIR': image.select('B4'),
            'green': image.select('B2')}
                                          )).rename('GNDVI')
        msavi2 = ee.Image(image.expression('(2 * NIR + 1 - ((2 * NIR + 1)**2 - 8 * (NIR -Red))**(1/2))/2', {
            'NIR': image.select('B4'),
            'Red': image.select('B3')
        })).rename('MSAVI2');
        albedo = image.expression('(0.356*B+0.130*R+0.373*NIR+0.085*SW1+0.072*SW2-0.018)',
                                  {'B': image.select('B1'), 'R': image.select('B3'),
                                   'NIR': image.select('B4'), 'SW1': image.select('B5'),
                                   'SW2': image.select('B7')}).rename('alfa')
        image = image.select(['B1','B2','B3','B4','B5','B7']).rename(['B2','B3','B4','B5','B6','B7'])

        return image.addBands(ndvi).addBands(nbr).addBands(nbr2).addBands(bai).addBands(gemi).addBands(mirbi).addBands(ndmi).addBands(savi).addBands(evi).addBands(gndvi).addBands(msavi2).addBands(albedo).copyProperties(img, ['system:time_start', 'system:time_end', 'LANDSAT_ID'])

#Cargar Fusion Table, base de datos incendios. Subir al servidor de gee
#table  = ee.FeatureCollection("users/italomoletto/MAPAS/Perimetros_IFM_CONAF_2015-2016");
#table  = ee.FeatureCollection("users/italomoletto/MAPAS/Perimetros_IFM_CONAF_2016-2017");
table = ee.FeatureCollection("users/italomoletto/MAPAS/FireDB_1985-2018_31052019_latlon")
#table = ee.FeatureCollection("users/italomoletto/MAPAS/fire_subset");
#table = table = ee.FeatureCollection("users/italomoletto/MAPAS/BD_incendios")
#Generar lista de incendios
lista = table.toList(table.size());
lista_for = len(lista.getInfo())
print(lista_for);

th_nbr = 0.1
th_ndvi1 = 0.2
th_ndvi2 = 0.2


L5_col = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA').map(maskClouds).map(get_INDEX_L57).map(adddate)
L7_col = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterDate('1999-01-01','2003-05-31').map(maskClouds).map(get_INDEX_L57).map(adddate)
L8_col = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').map(maskClouds).map(get_INDEX_L8).map(adddate)
L57_col = L5_col.merge(L7_col)
sat = L57_col.merge(L8_col).sort('system:time_start')


for i in range(0,lista_for):
                print(i)
                #Crear Feature de cada incendio para su tratamiento individual
                inc = ee.Feature(lista.get(i));     # seleccion del incendio a tratar(pendiente el 10,13)
                #print(inc);
            
                # rescatar las coordenadas de la tabla de atributos y generar el punto como un ee.Geometry.Point()
                inc2 = inc.centroid().geometry()
                #lat = inc.get('lat');
                #long = inc.get('long');
                #inc2 = ee.Geometry.Point([long, lat]); # creacion de punto central
                   
                #Umbral RdNBR (valor umbral del indice para ser considerado como area quemada)
                #   ind_RdNBR = 10;#10; # setear umbral de RdNBR
                
                # Generación de área de análisis en función al tamaño conocido del incendio. 
                distance = ee.Number(inc.get('sup_GIS'));
                distance = ee.Number(ee.Algorithms.If(distance, distance, ee.Number(1000)));  # condicion si el area estimada por conaf es cero
            
                dist2 = distance.log().multiply(2000); #multiply(20).add(1000); # operacion para indicar tamaño del buffer
                region = inc2.buffer(dist2); # buffer que contiene el area quemada             
              
                #valor del buffer generado para filtrar por distancia al polígono quemado de mayor superficie
                del_buf = distance.log().multiply(100); 
                #print('distancia buffer',del_buf);
            
              
              
                ##############/--- BOX 2 ---###########################/
                ################################################/
                ################################################/
                # Selección de colección de imagenes y generación de condición pre y post incendio#######
                ################################################/
                ################################################/
                
                #Seleccionar fecha para analisis
                #Fecha de inicio de incendio

                diai = inc.get('iniciodia');                                      
                mesi = inc.get('iniciomes');                                       
                anoi = inc.get('inicioano');                                       
            
                #Fecha de control de incendio 
                diaf = inc.get('controldia');                                   
                mesf = inc.get('controlmes');                                    
                anof = inc.get('controlano');
                print(diaf.getInfo())
                print(mesf.getInfo())
                print(anof.getInfo())
                #################################################
                #Fecha de inicio y termino del incendio en formato ee

                inicio = str(anoi.getInfo()) + '-' + str(mesi.getInfo()) + '-' + str(diai.getInfo())
                fin    = str(anof.getInfo()) + '-' + str(mesf.getInfo()) + '-' + str(diaf.getInfo())

                #inicio = ee.Date(inc.get('Inicio')) #Brull 2016
                #inicio = ee.Date(inc.get('inicio')) #Brull 2017
                #fin = inicio.advance(ee.Number(1),'month') #valor arbitrario
                inicio_th = ee.Date(inicio).advance(ee.Number(-1),'year');
                fin_th = ee.Date(fin)
            
                #Definir umbral de fechas pre y post incendio a ser evaluada, 86400000 = milisegundos por dia (24*60*60*1000)
                try:
                    inib = ee.Date(inicio).millis().getInfo() - 10000000000;
                    inibf = ee.Date(inib);
                    finb = ee.Date(fin).millis().getInfo() + 10000000000;
                    finbf = ee.Date(finb);
                except ee.EEException:
                    print('no valid date found in fire point')
                    continue
            
                ##############/--- BOX 3 ---###########################/ 
                ################################################/
                ################################################/
                # Generación de composición de imagenes e imagen final pre y post incendio ##########
                ################################################/
                ################################################/
            
                #Colección de imagenes Pre incendio TOA
                mosaicpre = ee.ImageCollection(sat).filterBounds(region).filterDate(inibf,inicio).map(maskClouds).map(adddate).sort('system:time_start', True)#.map(get_INDEX_L8);
               
                #Colección de imagenes Post incendio TOA
                mosaicpos = ee.ImageCollection(sat).filterBounds(region).filterDate(fin,finbf).map(maskClouds).map(adddate).sort('system:time_start', False)#.map(get_INDEX_L8);
            
                #Colección de imágenes año Post indencio y pre incendio para regla de decisión
                max_NDVI = ee.ImageCollection(sat).filterDate(inicio_th,fin_th).map(get_INDEX_L8).select('NDVI').max()
                min_NBR = ee.ImageCollection(sat).filterDate(inicio_th,fin_th).map(get_INDEX_L8).select('NBR').min()
            
                # crer mosaico (prioriza la imagen mas reciente, que seria la ultima del imagecollection, por lo tanto hay que invertir con sort para que use la mas antigua en la post)
                mosaic_pre = mosaicpre.mosaic();
                mosaic_pos = mosaicpos.mosaic();
            
                # Generación imagen pre incendio recortada al area de estudio
                pref = mosaic_pre.clip(region);
                postf = mosaic_pos.clip(region);

                ##############/--- BOX 4 ---###########################/ 
                ################################################/
                ################################################/
                # Generación de indices NBR pre y post incendio (dNBR) ####################/
                ################################################/
                ################################################/
                #Crear NDVI para filtrar cicatriz y definir el umbral (prueba)
                ndvipos = postf.select('NDVI');
                ndvipre = pref.select('NDVI');
            
            
                #AWEI
                aweipos = pref.expression('4*(GREEN-SWIR1)-(0.25*NIR+2.75*SWIR2)',{'GREEN':postf.select('B3'),'SWIR1':postf.select('B6'),'SWIR2':postf.select('B7'),'NIR':postf.select('B5')})
                #shadow mask on slope
                dataset = ee.Image('USGS/SRTMGL1_003');
                elevation = dataset.select('elevation');
                slope = ee.Terrain.slope(elevation)
                slope_bool = slope.lt(5)
                aweipos_bool = aweipos.lt(0).add(slope_bool).neq(2)
                aweipos_bool_ndvi = aweipos_bool.multiply(ndvipre.lt(0.2)).eq(0)
                #Map.addLayer(aweipos_bool_ndvi)
                # aweipos = postf.expression('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2',{'GREEN':postf.select('B3'),'SWIR1':postf.select('B6'),'SWIR2':postf.select('B7'),'NIR':postf.select('B5'),'BLUE':postf.select('B2')})
                #print(aweipos)
            
                # Calculo de indice NBR
                pref_nbr = pref.select('NBR')
                postf_nbr = postf.select('NBR')
            
                # calcular RdNBR con NBR normalizado
                # RdNBR = pref_nbr.subtract(postf_nbr); # RdNBR con NBR normalizados
                #print('RdNBR datos', RdNBR.reduceRegion(ee.Reducer.min(), region, 30),RdNBR.reduceRegion(ee.Reducer.max(), region, 30));
                RdNBR = pref_nbr.expression("(NBRpre-NBRpost)/((NBRpre/1000)**(1/2))",{'NBRpre':pref_nbr.select('NBR'),'NBRpost':postf_nbr.select('NBR')}).rename('RdNBR');
                dNBR = pref_nbr.expression("(NBRpre-NBRpost)",{'NBRpre':pref_nbr.select('NBR'),'NBRpost':postf_nbr.select('NBR')}).rename('dNBR');
                image_export = mosaic_pre.select(['B2','B3','B4','B6','B7']).addBands(mosaic_pos.select(['B2','B3','B4','B6','B7'])).addBands(dNBR).addBands(ndvipre).addBands(ndvipos)
                #minimal value 5
                percentile = 90
                ind_RdNBR_init = ee.Image(ee.Number(RdNBR.reduceRegion(ee.Reducer.percentile([percentile]),region,30).get('RdNBR')))
                ind_RdNBR_cond = ind_RdNBR_init.gt(5) #Minimal empirical value
                ind_RdNBR_cond_remap = ind_RdNBR_cond.remap([1,0],[0,5])
                ind_RdNBR = ind_RdNBR_init.multiply(ind_RdNBR_cond).add(ind_RdNBR_cond_remap)
                #Map.addLayer(ind_RdNBR)
                #print(ind_RdNBR_init)
                #print(ind_RdNBR)
                threshold_ndvi1 = max_NDVI.gt(th_ndvi1)
                threshold_ndvi2 = max_NDVI.subtract(0.2).gt(th_ndvi2)
                threshold_nbr = min_NBR.subtract(postf_nbr).gt(th_nbr)
                contaminatedmask_1 = aweipos_bool_ndvi.And(threshold_ndvi1).And(threshold_ndvi2)
                contaminatedmask_2 =  aweipos_bool_ndvi.And(threshold_ndvi1).And(threshold_ndvi2).And(threshold_nbr)
                #Map.addLayer(aweipos_bool)
                #Map.addLayer(threshold_nbr)
                #Map.addLayer(contaminatedmask)
                firemask_and_contmask = RdNBR.gt(ind_RdNBR).And(contaminatedmask_1);#parametro iable en función a la latitud~precipitación
                firemask_and_contmask_2 = RdNBR.gt(ind_RdNBR).And(contaminatedmask_2);#parametro iable en función a la latitud~precipitación
                firemask = firemask_and_contmask.updateMask(firemask_and_contmask.eq(1));
                firemask_2 = firemask_and_contmask.updateMask(firemask_and_contmask_2.eq(1));
            
            
                ##############/--- BOX 5 ---###########################/ 
                ################################################/
                ################################################/
                # Creación de poligono y depuración de cicatriz de incendios##################
                ################################################/
                ################################################/
            
            
                #Transformar a vector 
                firevect = firemask.addBands(ndvipos).reduceToVectors( geometry= region, crs= ndvipos.projection(),scale= 30,geometryType= 'polygon', eightConnected= False,labelProperty= 'label',reducer= ee.Reducer.mean())
            
                # Incluir el conteo de pixeles por poligono en las propiedades de "feature collection"
                #calculo de superficie por pixel
                count = firemask.multiply(ee.Image.pixelArea());
              
                firevect_area = count.reduceRegions(collection=firevect,reducer= ee.Reducer.sum(),scale= 30)
                    #Map.addLayer(firevect)
            
                #Filtrar por superficie
                #Eliminar todos los poligonos con superficie menor a sup quemada / 50 en numero de pixels
                sup = distance.multiply(10000).divide(5);
            
                firevect_poly_area = firevect_area.filterMetadata('sum', 'greater_than', 10000)
                #Filtrar por respuesta espectral banda 5
                b5 = postf.select('B5');
              
                firevect_poly_area_b5 = b5.reduceRegions(collection=firevect_poly_area,reducer= ee.Reducer.min(),scale= 30)
                
                firevect_poly_area_b5f = firevect_poly_area_b5.filterMetadata('min', 'less_than', 0.15).filterMetadata('min', 'greater_than', 0.01);
            
                # Filtrar por ubicación
                # encontrar mayor superficie cerca del centroide
                bufsel0 = inc2.buffer(1000);  
                firevect_poly_area_b5f_buf0 = firevect_poly_area_b5f.filterBounds(bufsel0);
                #     firevect_poly_area_b5f_buf0 = firevect_area.filterBounds(bufsel0)
                #    Map.addLayer(firevect_ndvi_area_b5f_aguaf_buf0)
                # Filtrar por distncia al incendios mas grnade cercano al centro
                maxArea = firevect_poly_area_b5f_buf0.sort('sum', False).limit(1).geometry(); # obtener poligono de mayor superficie
              
                bufsel = maxArea.buffer(del_buf);  
                firevect_poly_area_b5f_buf = firevect_poly_area_b5f.filterBounds(bufsel);
                # firevect_poly_area_b5f_buf = firevect_poly_area_b5f_buf0.filterBounds(bufsel);
                # Calcular superficie estimada incendio
                area2 = firevect_poly_area_b5f_buf.reduceColumns(reducer= ee.Reducer.sum(),selectors= ['sum'])
            
                #CONDICIONAL PARA CONTROL DE CALIDAD DE LA IMAGEN EN CASO DE QUE SEA SUPERIOR A 2000 ha        
                fire_scar_img = ee.Image(ee.Algorithms.If(ee.Number(area2.get('sum')).gt(20000000),firemask_2,firemask))
                #Map.addLayer(fire_scar_img)
                #print(fire_scar_img)
                #PASO POS CONDICIONAL, REPETIR TRANSFORMACION A VECTOR
                
                #Transformar a vector 
                firevect = fire_scar_img.addBands(ndvipos).reduceToVectors(geometry= region,crs= ndvipos.projection(),scale= 30,geometryType= 'polygon',eightConnected= False,labelProperty= 'label',reducer= ee.Reducer.mean())
            
                # Incluir el conteo de pixeles por poligono en las propiedades de "feature collection"
                #calculo de superficie por pixel
                count = fire_scar_img.multiply(ee.Image.pixelArea());
              
                firevect_area = count.reduceRegions(collection=firevect,reducer= ee.Reducer.sum(),scale= 30)
                #Map.addLayer(firevect)
            
                #Filtrar por superficie
                #Eliminar todos los poligonos con superficie menor a sup quemada / 50 en numero de pixels
                sup = distance.multiply(10000).divide(5);
            
                firevect_poly_area = firevect_area.filterMetadata('sum', 'greater_than', 10000)
                #Filtrar por respuesta espectral banda 5
                b5 = postf.select('B5');
              
                firevect_poly_area_b5 = b5.reduceRegions(collection=firevect_poly_area,reducer= ee.Reducer.min(),scale= 30);
                
                firevect_poly_area_b5f = firevect_poly_area_b5.filterMetadata('min', 'less_than', 0.15).filterMetadata('min', 'greater_than', 0.01);
            
                # Filtrar por ubicación
                # encontrar mayor superficie cerca del centroide
                bufsel0 = inc2.buffer(1000);  
                firevect_poly_area_b5f_buf0 = firevect_poly_area_b5f.filterBounds(bufsel0);
                #     firevect_poly_area_b5f_buf0 = firevect_area.filterBounds(bufsel0)
                #    Map.addLayer(firevect_ndvi_area_b5f_aguaf_buf0)
                # Filtrar por distncia al incendios mas grnade cercano al centro
                maxArea = firevect_poly_area_b5f_buf0.sort('sum', False).limit(1).geometry(); # obtener poligono de mayor superficie
              
                bufsel = maxArea.buffer(del_buf);  
                firevect_poly_area_b5f_buf = firevect_poly_area_b5f.filterBounds(bufsel);
                # firevect_poly_area_b5f_buf = firevect_poly_area_b5f_buf0.filterBounds(bufsel);
                # Calcular superficie estimada incendio
                area2 = firevect_poly_area_b5f_buf.reduceColumns(reducer= ee.Reducer.sum(),selectors= ['sum']);
            
                #FILTRO POR KERNEL y PRESENCIA DE AGUA + PIXELES CONTAMINADOS
                """
                fire_scar_img_filtered = firevect_poly_area_b5f_buf
                .filter(ee.Filter.notNull(['label']))
                .reduceToImage({
                        properties: ['label'],
                        reducer: ee.Reducer.first()
                        }).updateMask(contaminatedmask);
                    
                    
                fire_scar_img_filtered_kernel = fire_scar_img_filtered.unmask(0).reduceNeighborhood({
                        reducer: ee.Reducer.mean(),
                        kernel: ee.Kernel.square(30,'meters'),
                        }).gt(0.9).remap([1],[1]);
                    
                #Transformar a vector 
                firevect_filtered = fire_scar_img_filtered.addBands(ndvipos).reduceToVectors({ 
                        geometry: region,    
                        crs: ndvipos.projection(),
                        scale: 30,
                        geometryType: 'polygon',
                        eightConnected: false,
                        labelProperty: 'label',
                        reducer: ee.Reducer.mean(),
                });
            
            
                # fire_scar_vect_filtered_kernel_area = fire_scar_img_filtered_kernel.multiply(ee.Image.pixelArea());
                fire_scar_img_filtered_area = fire_scar_img_filtered.multiply(ee.Image.pixelArea());
                #fire_scar_vect_filtered_kernel_area
                fire_scar_vect_filtered = fire_scar_img_filtered_area.reduceRegions({
                collection:firevect_filtered,
                reducer: ee.Reducer.sum(),
                scale: 30,
                });
                """
                # fire_scar_vect_filtered_kernel_filt2 = fire_scar_vect_filtered_kernel.filterMetadata('sum', 'greater_than', 10000)
                fire_scar_vect_filtered_filt2 = firevect_poly_area_b5f_buf.filterMetadata('sum', 'greater_than', 10000)
            
                mean_RdNBR = RdNBR.reduceRegions(collection=fire_scar_vect_filtered_filt2,reducer= ee.Reducer.mean(),scale= 30) 
            
                meanIntense = mean_RdNBR.sort('mean', False).limit(1).geometry(); # obtener poligono de mayor superficie
                bufsel_rdnbr = meanIntense.buffer(150); #Buffer de 5 pixeles 
                fire_scar_vect_filtered_filt2_areabuf =fire_scar_vect_filtered_filt2.filterBounds(bufsel_rdnbr)
                #print('intensidad buffer',mean_RdNBR.sort('mean', false))
                #Map.addLayer(fire_scar_vect_filtered_filt2)
                #Map.addLayer(bufsel_rdnbr)
            
                area3 = fire_scar_vect_filtered_filt2_areabuf.reduceColumns(reducer= ee.Reducer.sum(),selectors= ['sum'])
                #print('area incendiada', area3)
                #fire_scar_vect_filtered_filt2
                #print('area',area2);
                #print(mosaicpre,mosaicpos)
                ##############/--- BOX 7 ---###########################/ 
                ################################################/
                ################################################/
                # Exportación de datos ####################################/
                ################################################/
                ################################################/ 
                
                #Exportar shape
                scar = inc.get('KEY').getInfo();
                #scar = inc.get('Incendio') #brull 2016
                #scar = inc.get('incendio') #brull 2017
                #scars = str(scar.getInfo()).replace('Ñ','N').replace('//','').replace('/','')
                #scars = str(int(scar.getInfo()));
                id = "Fire_scars_ID"; 

                      #print(dwl);
                print('exporting data '+id + str(scar) + ' Date ' + fin+' to folder: '+fin[0:4])
                #fecha = datetime.datetime.utcfromtimestamp(int(inicio.getInfo()['value']) / 1000).strftime('%Y-%m-%d')
                #print(fecha)
                export = ee.batch.Export.image.toDrive(\
                                                       image=image_export,
                                                       description= id + str(scar) + '_Date_' + fin, #+ '_u' + ind_RdNBR,
                                                       maxPixels= 7865781687,
                                                       region=region.getInfo()['coordinates'],
                                                       folder = 'FireScar-CL_'+fin[0:4],
                                                       scale = 30,
                                                       fileFormat='GeoTIFF')#,
                export.start()
                
                #fire_scar_vect_filtered_filt2_areabuf

                export = ee.batch.Export.table.toDrive(\
                                                       collection=firevect_poly_area_b5f_buf,
                                                       description= id+'_' + str(scar) + '_Date_' + fin,# + '_u' + ind_RdNBR,
                                                       folder= 'FireScar-CL_'+fin[0:4],
                                                       fileFormat= 'KML')
                export.start()