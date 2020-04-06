//--------------------------------------------------------------------------------------------------------------------------------
//--Este script carga la base de datos de incendios de donde se obtiene la fecha de inicio y termino del incendio, ---------------
//--su ubicación estimada y su superficie estimada. Selecciona las imagenes pre y post incendio mediante una reducción de media,--
//--Crea un área de analisis en donde compara el pre y post indice NBR.-----------------------------------------------------------
//--Genera el perimetro de cada incendio y su severidad.-------------------------------------------------------------------------- 
//--------------------------------------------------------------------------------------------------------------------------------


/////////////////////////////--- BOX 0 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Definir funciones para filtrar las imagenes a utiliziar - Seleccion de datos satelitales //////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Generar mascara de nubes y sombra de nubes
// Obtención del valor de los bits correspondiente a los valores por pixel de la calidad de las imagenes
                        var getQABits = function(image, start, end, newName) {
                            // Calculo de bits a extraer
                            var pattern = 0;
                            for (var i = start; i <= end; i++) {
                               pattern += Math.pow(2, i);
                            }
                            // Generación de banda de calidad por pixel en bit
                            // nuevo nombre a imagen
                            return image.select([0], [newName])
                                          .bitwiseAnd(pattern)
                                          .rightShift(start);
                            };

// Función para enmascarar sombra de nubes y nieve (7,8 sombra de nube y 9,10 nieve y hielo)
                          var cloud_shadows = function(image) {
                            // Seleccionar la banda de calidad
                            var QA = image.select(['BQA']);
                            // Aplicar función 
                            return getQABits(QA, 7,8, 'Cloud_shadows').eq(1);
                            };
                          
                          // Función para enmascarar nubes
                            var clouds = function(image) {
                            // Seleccionar la banda de calidad
                            var QA = image.select(['BQA']);
                            // Aplicar función 
                            return getQABits(QA, 4,4, 'Cloud').eq(0);
                            };
                            
                            // Función para enmascarar nieve
                          var snow = function(image) {
                            // Seleccionar la banda de calidad
                          var QA = image.select(['BQA']);
                          // Aplicar función 
                           return getQABits(QA, 9,10, 'Snow').eq(1);
                           };
                          
                          // Generación función de mascara
                          var maskClouds = function(image) {
                            var cs = cloud_shadows(image);
                            var c = clouds(image);
                            var s = snow(image);
                            image = image.updateMask(cs);
                            image = image.updateMask(s);
                            return image.updateMask(c);
                            
                            };
                            
// Agregar banda representando la fecha de adquisicion //https://developers.google.com/earth-engine/ic_composite_mosaic
                            var adddate = function(image) {
                              return image.addBands(image.metadata('system:time_start'));
                            };

var get_INDEX_L8 = function(img) {
    var image = img.divide(10000)
    var ndvi = ee.Image(image.expression('(NIR - Red) / (NIR + Red)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4')}
                                     )).rename('NDVI')
    var nbr = ee.Image(image.expression('(NIR - SWIR2) / (NIR + SWIR2)', {
        'NIR': image.select('B5'),
        'SWIR2': image.select('B7')}
                                     )).rename('NBR')
    var savi = ee.Image(image.expression('((NIR - Red) / (NIR + Red + L)) * (1 + L)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4'),
        'L': 0.5, }
                                     )).rename('SAVI')

    var evi = ee.Image(image.expression('2.5 * (NIR - Red) / ((NIR + C1 * Red - C2 * Blue) + L)', {
        'NIR': image.select('B5'),
        'Red': image.select('B4'),
        'C1': 6,
        'C2': 7.5,
        'Blue': image.select('B2'),
        'L': 1}
                                    )).rename('EVI')

    var gndvi = ee.Image(image.expression('(NIR - green) / (NIR + green)', {
        'NIR': image.select('B5'),
        'green': image.select('B3')}
                                      )).rename('GNDVI')

    var msavi2 = ee.Image(image.expression('(2 * NIR + 1 - ((2 * NIR + 1)**2 - 8 * (NIR -Red))**(1/2))/2', {
        'NIR': image.select('B5'),
        'Red': image.select('B4')
    })).rename('MSAVI2');

    var albedo = image.expression('(0.356*B+0.130*R+0.373*NIR+0.085*SW1+0.072*SW2-0.018)',
                              {'B': image.select('B2'), 'R': image.select('B4'),
                               'NIR': image.select('B5'), 'SW1': image.select('B6'),
                               'SW2': image.select('B7')}).rename('alfa')

    return image.addBands(ndvi).addBands(nbr).addBands(savi).addBands(evi).addBands(gndvi).addBands(msavi2).addBands(albedo).copyProperties(img, ['system:time_start', 'system:time_end','LANDSAT_ID'])
}
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //Selección de colección de imagenes a utilizar en función del año del incendio
  //Para el año 2013 o más reciente utilizar Landsat 8 OLI
    var satellite = 'Landsat 8 OLI';
    var sat = 'LANDSAT/LC08/C01/T1_TOA';
    
  // Definir bandas de acuerdo al sensor utilizado B5 = NIR, B7 = SWIR2
    var B4 = 'B5';
    var B7 = 'B7';
  
/////////////////////////////--- BOX 1 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Importar incendios, creación de tabla y definición área de analisis //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//Cargar Fusion Table, base de datos incendios. Subir al servidor de gee
//var table = ee.FeatureCollection("users/italomoletto/MAPAS/fire_subset");
var table = ee.FeatureCollection("users/italomoletto/MAPAS/BD_incendios")
//Generar lista de incendios
  var lista = table.toList(table.size());
  print(lista);

//##@@ --------------------------------------Comienzo del codigo que entra al loop ----------------------------------------@@##//

//Crear Feature de cada incendio para su tratamiento individual
  var inc = ee.Feature(lista.get(944));     // seleccion del incendio a tratar(pendiente el 10,13)
  var inc = table.filterMetadata('KEY','equals',189210).first() 
  print(inc);

// rescatar las coordenadas de la tabla de atributos y generar el punto como un ee.Geometry.Point()
    var lat = inc.get('lat');
    var long = inc.get('long');
  var inc2 = ee.Geometry.Point([long, lat]); // creacion de punto central
       
//Umbral RdNBR (valor umbral del indice para ser considerado como area quemada)
//  var ind_RdNBR = 10;//10; // setear umbral de RdNBR

// Generación de área de análisis en función al tamaño conocido del incendio. 
    var distance = ee.Number(inc.get('area'));
    var distance = ee.Number(ee.Algorithms.If(distance, distance, ee.Number(1000)));  // condicion si el area estimada por conaf es cero

  var dist2 = distance.log().multiply(2000); //multiply(20).add(1000); // operacion para indicar tamaño del buffer
  var region = inc2.buffer(dist2); // buffer que contiene el area quemada             
  
//valor del buffer generado para filtrar por distancia al polígono quemado de mayor superficie
  var del_buf = distance.log().multiply(100); 
  print('distancia buffer',del_buf);

  
  
/////////////////////////////--- BOX 2 ---///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Selección de colección de imagenes y generación de condición pre y post incendio//////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//Seleccionar fecha para analisis
//Fecha de inicio de incendio
  var diai = inc.get('iniciodia');                                      
  var mesi = inc.get('iniciomes');                                       
  var anoi = inc.get('inicioano');                                       

//Fecha de control de incendio 
  var diaf = inc.get('controldia');                                   
  var mesf = inc.get('controlmes');                                    
  var anof = inc.get('controlano');                                     

//////////////////////////////////////////////////////////////////////////////////////////////////
//Fecha de inicio y termino del incendio en formato ee 
  var inicio = anoi.getInfo() + '-' + mesi.getInfo() + '-' + diai.getInfo();
  var fin    = anof.getInfo() + '-' + mesf.getInfo() + '-' + diaf.getInfo();
  
  var inicio_th = ee.Date(inicio).advance(ee.Number(-1),'year');
  var fin_th = ee.Date(fin)
print(inicio_th)
print(fin_th)
//Definir umbral de fechas pre y post incendio a ser evaluada, 86400000 = milisegundos por dia (24*60*60*1000)
    var inib = ee.Date(inicio).millis().getInfo() - 10000000000; 
  var inibf = ee.Date(inib);
    var finb = ee.Date(fin).millis().getInfo() + 10000000000;
  var finbf = ee.Date(finb);                     

/////////////////////////////--- BOX 3 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Generación de composición de imagenes e imagen final pre y post incendio ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//Colección de imagenes Pre incendio TOA
  var mosaicpre = ee.ImageCollection(sat)
  .filterBounds(region)
  .filterDate(inibf,inicio)
  .map(maskClouds)
  .map(adddate)
  .sort('system:time_start', true)//.map(get_INDEX_L8);
   
//Colección de imagenes Post incendio TOA
  var mosaicpos = ee.ImageCollection(sat)
  .filterBounds(region)
  .filterDate(fin,finbf)
  .map(maskClouds)
  .map(adddate)
  .sort('system:time_start', false)//.map(get_INDEX_L8);

//Colección de imágenes año Post indencio y pre incendio para regla de decisión
var max_NDVI = ee.ImageCollection(sat).filterDate(inicio_th,fin_th).map(get_INDEX_L8).select('NDVI').max()
var min_NBR = ee.ImageCollection(sat).filterDate(inicio_th,fin_th).map(get_INDEX_L8).select('NBR').min()

// crer mosaico (prioriza la imagen mas reciente, que seria la ultima del imagecollection, por lo tanto hay que invertir con sort para que use la mas antigua en la post)
  var mosaic_pre = mosaicpre.mosaic();
  var mosaic_pos = mosaicpos.mosaic();

// Generación imagen pre incendio recortada al area de estudio
  var pref = mosaic_pre.clip(region);
  var postf = mosaic_pos.clip(region);
  
/////////////////////////////--- BOX 4 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Generación de indices NBR pre y post incendio (dNBR) /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//Crear NDVI para filtrar cicatriz y definir el umbral (prueba)
  var ndvipos = postf.normalizedDifference(['B5','B4']);
  var ndvipre = pref.normalizedDifference(['B5','B4']);


//AWEI
  var aweipos = pref.expression('4*(GREEN-SWIR1)-(0.25*NIR+2.75*SWIR2)',{'GREEN':postf.select('B3'),'SWIR1':postf.select('B6'),'SWIR2':postf.select('B7'),'NIR':postf.select('B5')})
//shadow mask on slope
var dataset = ee.Image('USGS/SRTMGL1_003');
var elevation = dataset.select('elevation');
var slope = ee.Terrain.slope(elevation)
var slope_bool = slope.lt(5)
var aweipos_bool = aweipos.lt(0).add(slope_bool).neq(2)
var aweipos_bool_ndvi = aweipos_bool.multiply(ndvipre.lt(0.2)).eq(0)
//Map.addLayer(aweipos_bool_ndvi)
  //var aweipos = postf.expression('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2',{'GREEN':postf.select('B3'),'SWIR1':postf.select('B6'),'SWIR2':postf.select('B7'),'NIR':postf.select('B5'),'BLUE':postf.select('B2')})
  //print(aweipos)

// Calculo de indice NBR
  var pref_nbr = pref.expression('float(b("'+B4+'") - b("'+B7+'")) / (b("'+B4+'") + b("'+B7+'"))');
  var postf_nbr = postf.expression('float(b("'+B4+'") - b("'+B7+'")) / (b("'+B4+'") + b("'+B7+'"))');

// calcular RdNBR con NBR normalizado
  //var RdNBR = pref_nbr.subtract(postf_nbr); // RdNBR con NBR normalizados
    //print('RdNBR datos', RdNBR.reduceRegion(ee.Reducer.min(), region, 30),RdNBR.reduceRegion(ee.Reducer.max(), region, 30));
  var firesev = pref_nbr.addBands(postf_nbr);
   var RdNBR = firesev.expression("(b('B5') - b('B5_1'))/(sqrt(abs(b('B5')/1000)))").toFloat().rename('RdNBR'); // metodo anterior
  //minimal value 5
  var percentile = 90
  var RdNBR_threshold = ee.Number(5)
  var ind_RdNBR_init = ee.Image(ee.Number(RdNBR.reduceRegion(ee.Reducer.percentile([percentile]),region,30).get('RdNBR')))
  var ind_RdNBR_cond = ind_RdNBR_init.gt(RdNBR_threshold) //Minimal empirical value
  var ind_RdNBR_cond_remap = ind_RdNBR_cond.remap([1,0],[0,RdNBR_threshold])
  var ind_RdNBR = ind_RdNBR_init.multiply(ind_RdNBR_cond).add(ind_RdNBR_cond_remap)
  //Map.addLayer(ind_RdNBR)
  print(ind_RdNBR_init)
  print(ind_RdNBR)
  var ind_RdNBR_number =  ee.Number(ind_RdNBR.reduceRegion(ee.Reducer.mean(),region,30).get('constant'))
  print(ind_RdNBR_number,'mean value of rdnbr')
  var threshold_ndvi1 = max_NDVI.gt(0.2)
  var threshold_ndvi2 = max_NDVI.subtract(0.2).gt(0.2)
  var threshold_nbr = min_NBR.subtract(postf_nbr).gt(0.1)
  var contaminatedmask_1 = aweipos_bool_ndvi.and(threshold_ndvi1).and(threshold_ndvi2)
  var contaminatedmask_2 =  aweipos_bool_ndvi.and(threshold_ndvi1).and(threshold_ndvi2).and(threshold_nbr)
  //Map.addLayer(aweipos_bool)
  Map.addLayer(threshold_nbr)
  //Map.addLayer(contaminatedmask_2)
  //var firemask_and_contmask = RdNBR.gt(ind_RdNBR).and(contaminatedmask_1);//parametro variable en función a la latitud~precipitación
  //var firemask_and_contmask_2 = RdNBR.gt(ind_RdNBR).and(contaminatedmask_2);//parametro variable en función a la latitud~precipitación
  var firemask_and_contmask =  ee.Image(ee.Algorithms.If(ind_RdNBR_number.gt(12),RdNBR.gt(RdNBR_threshold),RdNBR.gt(ind_RdNBR).and(contaminatedmask_1))).add(RdNBR.gt(20)).gte(1)
  var firemask_and_contmask_2 =  ee.Image(ee.Algorithms.If(ind_RdNBR_number.gt(12),RdNBR.gt(RdNBR_threshold),RdNBR.gt(ind_RdNBR).and(contaminatedmask_2))).add(RdNBR.gt(20)).gte(1)
  Map.addLayer(firemask_and_contmask_2)
  var firemask = firemask_and_contmask.updateMask(firemask_and_contmask.eq(1));
  var firemask_2 = firemask_and_contmask.updateMask(firemask_and_contmask_2.eq(1));
  //Map.addLayer(firemask_2)

/////////////////////////////--- BOX 5 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Creación de poligono y depuración de cicatriz de incendios////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


//Transformar a vector 
  var firevect = firemask.addBands(ndvipos).reduceToVectors({ 
   geometry: region,    
   crs: ndvipos.projection(),
   scale: 30,
   geometryType: 'polygon',
   eightConnected: false,
   labelProperty: 'label',
   reducer: ee.Reducer.mean(),
   });

// Incluir el conteo de pixeles por poligono en las propiedades de "feature collection"
  //calculo de superficie por pixel
  var count = firemask.multiply(ee.Image.pixelArea());
  
  var firevect_area = count.reduceRegions({
    collection:firevect,
    reducer: ee.Reducer.sum(),
    scale: 30,
    });
//Map.addLayer(firevect)

//Filtrar por superficie
//Eliminar todos los poligonos con superficie menor a sup quemada / 50 en numero de pixels
  var sup = distance.multiply(10000).divide(5);

  var firevect_poly_area = firevect_area.filterMetadata('sum', 'greater_than', 10000)
//Filtrar por respuesta espectral banda 5
  var b5 = postf.select('B5');
  
  var firevect_poly_area_b5 = b5.reduceRegions({
                              collection:firevect_poly_area,
                              reducer: ee.Reducer.min(),
                              scale: 30,
                              });
    
  var firevect_poly_area_b5f = firevect_poly_area_b5.filterMetadata('min', 'less_than', 0.15)
                               .filterMetadata('min', 'greater_than', 0.01);

// Filtrar por ubicación
  // encontrar mayor superficie cerca del centroide
    var bufsel0 = inc2.buffer(1000);  
    var firevect_poly_area_b5f_buf0 = firevect_poly_area_b5f.filterBounds(bufsel0);
//    var firevect_poly_area_b5f_buf0 = firevect_area.filterBounds(bufsel0)
//    Map.addLayer(firevect_ndvi_area_b5f_aguaf_buf0)
  // Filtrar por distncia al incendios mas grnade cercano al centro
      var maxArea = firevect_poly_area_b5f_buf0.sort('sum', false).limit(1).geometry(); // obtener poligono de mayor superficie
    var bufsel = maxArea.buffer(del_buf);  
    var firevect_poly_area_b5f_buf = firevect_poly_area_b5f.filterBounds(bufsel);
    //var firevect_poly_area_b5f_buf = firevect_poly_area_b5f_buf0.filterBounds(bufsel);
// Calcular superficie estimada incendio
  var area2 = firevect_poly_area_b5f_buf.reduceColumns({
        reducer: ee.Reducer.sum(),
        selectors: ['sum'],
        });

//CONDICIONAL PARA CONTROL DE CALIDAD DE LA IMAGEN EN CASO DE QUE SEA SUPERIOR A 2000 ha        
var fire_scar_img = ee.Image(ee.Algorithms.If(ee.Number(area2.get('sum')).gt(20000000),firemask_2,firemask))
//Map.addLayer(fire_scar_img)
print(ee.Number(area2.get('sum')),'suma_condicional_area')
print(fire_scar_img,'area_incendiada')
//PASO POS CONDICIONAL, REPETIR TRANSFORMACION A VECTOR

//Transformar a vector 
  var firevect = fire_scar_img.addBands(ndvipos).reduceToVectors({ 
   geometry: region,    
   crs: ndvipos.projection(),
   scale: 30,
   geometryType: 'polygon',
   eightConnected: false,
   labelProperty: 'label',
   reducer: ee.Reducer.mean(),
   });

// Incluir el conteo de pixeles por poligono en las propiedades de "feature collection"
  //calculo de superficie por pixel
  var count = fire_scar_img.multiply(ee.Image.pixelArea());
  
  var firevect_area = count.reduceRegions({
    collection:firevect,
    reducer: ee.Reducer.sum(),
    scale: 30,
    });
//Map.addLayer(firevect)

//Filtrar por superficie
//Eliminar todos los poligonos con superficie menor a sup quemada / 50 en numero de pixels
  var sup = distance.multiply(10000).divide(5);

  var firevect_poly_area = firevect_area.filterMetadata('sum', 'greater_than', 10000)
//Filtrar por respuesta espectral banda 5
  var b5 = postf.select('B5');
  
  var firevect_poly_area_b5 = b5.reduceRegions({
                              collection:firevect_poly_area,
                              reducer: ee.Reducer.min(),
                              scale: 30,
                              });
    
  var firevect_poly_area_b5f = firevect_poly_area_b5.filterMetadata('min', 'less_than', 0.15)
                               .filterMetadata('min', 'greater_than', 0.01);

// Filtrar por ubicación
  // encontrar mayor superficie cerca del centroide
    var bufsel0 = inc2.buffer(1000);  
    var firevect_poly_area_b5f_buf0 = firevect_poly_area_b5f.filterBounds(bufsel0);
//    var firevect_poly_area_b5f_buf0 = firevect_area.filterBounds(bufsel0)
//    Map.addLayer(firevect_ndvi_area_b5f_aguaf_buf0)
  // Filtrar por distncia al incendios mas grnade cercano al centro
      var maxArea = firevect_poly_area_b5f_buf0.sort('sum', false).limit(1).geometry(); // obtener poligono de mayor superficie
  
    var bufsel = maxArea.buffer(del_buf);  
    var firevect_poly_area_b5f_buf = firevect_poly_area_b5f.filterBounds(bufsel);
    //var firevect_poly_area_b5f_buf = firevect_poly_area_b5f_buf0.filterBounds(bufsel);
print(firevect_poly_area_b5f.reduceColumns({
        reducer: ee.Reducer.sum(),
        selectors: ['sum'],
        }));
// Calcular superficie estimada incendio
  var area2 = firevect_poly_area_b5f_buf.reduceColumns({
        reducer: ee.Reducer.sum(),
        selectors: ['sum'],
        });

//FILTRO POR KERNEL y PRESENCIA DE AGUA + PIXELES CONTAMINADOS
/*
var fire_scar_img_filtered = firevect_poly_area_b5f_buf
  .filter(ee.Filter.notNull(['label']))
  .reduceToImage({
    properties: ['label'],
    reducer: ee.Reducer.first()
}).updateMask(contaminatedmask);


var fire_scar_img_filtered_kernel = fire_scar_img_filtered.unmask(0).reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.square(30,'meters'),
}).gt(0.9).remap([1],[1]);

//Transformar a vector 
  var firevect_filtered = fire_scar_img_filtered.addBands(ndvipos).reduceToVectors({ 
   geometry: region,    
   crs: ndvipos.projection(),
   scale: 30,
   geometryType: 'polygon',
   eightConnected: false,
   labelProperty: 'label',
   reducer: ee.Reducer.mean(),
   });


//var fire_scar_vect_filtered_kernel_area = fire_scar_img_filtered_kernel.multiply(ee.Image.pixelArea());
var fire_scar_img_filtered_area = fire_scar_img_filtered.multiply(ee.Image.pixelArea());
//fire_scar_vect_filtered_kernel_area
var fire_scar_vect_filtered = fire_scar_img_filtered_area.reduceRegions({
    collection:firevect_filtered,
    reducer: ee.Reducer.sum(),
    scale: 30,
    });
*/
//var fire_scar_vect_filtered_kernel_filt2 = fire_scar_vect_filtered_kernel.filterMetadata('sum', 'greater_than', 10000)
var fire_scar_vect_filtered_filt2 = firevect_poly_area_b5f_buf.filterMetadata('sum', 'greater_than', 10000)

var mean_RdNBR = RdNBR.reduceRegions({
                              collection:fire_scar_vect_filtered_filt2,
                              reducer: ee.Reducer.mean(),
                              scale: 30,
                              }); 

var meanIntense = mean_RdNBR.sort('mean', false).limit(1).geometry(); // obtener poligono de mayor superficie
var bufsel_rdnbr = meanIntense.buffer(150); //Buffer de 5 pixeles 
var fire_scar_vect_filtered_filt2_areabuf =fire_scar_vect_filtered_filt2.filterBounds(bufsel_rdnbr)
print('intensidad buffer',mean_RdNBR.sort('mean', false))
//Map.addLayer(fire_scar_vect_filtered_filt2)
//Map.addLayer(bufsel_rdnbr)

var area3 = fire_scar_vect_filtered_filt2_areabuf.reduceColumns({
        reducer: ee.Reducer.sum(),
        selectors: ['sum'],
        });
print('area incendiada', area3)
//fire_scar_vect_filtered_filt2
print('area',area2);
//print(mosaicpre,mosaicpos)
/////////////////////////////--- BOX 7 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Exportación de datos /////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////// 

//Exportar shape
  var scar = inc.get('KEY');
  var scars = scar.getInfo();
  var id = "Fire_scars_ID"; 

var dwl = RdNBR.addBands(firesev)
          .addBands(ndvipre)
          .addBands(ndvipos);
//print(dwl);

  Export.image.toDrive({
    image: dwl, //RdNBR
    description: id + scars + '_Date_' + inicio, //+ '_u' + ind_RdNBR,
    folder: 'Google engine',
    scale: 30,
    region: region,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });

  Export.table.toDrive({
    collection: fire_scar_vect_filtered_filt2_areabuf,
    description: id + scars + '_Date_' + inicio,// + '_u' + ind_RdNBR,
    folder: 'Google engine',
    fileFormat: 'KML'
  });
 
 //##@@ ------------------------------------------- Fin del codigo que entra al loop ---------------------------------------------@@##//
  
/////////////////////////////--- BOX 8 ---/////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// Desplegar mapas /////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
  
//Centrar área de trabajo en la ubicación reportada del incendio
  Map.centerObject(inc,13);   
  
// generar imagen pre y post fuego
  Map.addLayer(pref, {bands: ['B6', 'B5', 'B4'],min: [0,0,0],max: [0.35, 0.35, 0.35]}, 'img pref',false); 
  Map.addLayer(postf, {bands: ['B6', 'B5', 'B4'],min: [0,0,0],max: [0.35, 0.35, 0.35]}, 'img post',false);
  
//Generar imagen de RdNBR
  //var dnbrPar = {min: -1, max: 2, palette: ['darkgreen','green','yellow','orange','red','darkred']};
  //Map.addLayer(pref_nbr,dnbrPar,'NBR prefire', false);
  //Map.addLayer(postf_nbr,dnbrPar,'NBR postfire', false);
  Map.addLayer(RdNBR, {min: 0, max: 25, palette: ['darkgreen','green','yellow','orange','red','darkred']}, 'RdNBR');
  //Map.addLayer(dnbr,{min: 0, max: 1, palette: ['darkgreen','green','yellow','orange','red','darkred']},'dNBR norm');

print(fire_scar_vect_filtered_filt2)
print(fire_scar_vect_filtered_filt2_areabuf)
// Mascaras de zonas quemadas en base a umbral definido
  //Map.addLayer(firemask, {}, 'Fire Mask', false);
  Map.addLayer(firevect_poly_area_b5f_buf,{},'huella final shp',false);
  //Map.addLayer(fire_scar_vect_filtered_kernel_filt2,{},'huella fixed kernel')
  Map.addLayer(fire_scar_vect_filtered_filt2_areabuf,{},'scar qa fixed')
//Desplegar puntos en el mapa
  Map.addLayer(table, {}, 'Puntos Incendios', false);

// Desplegar filtro por ndvi
  //Map.addLayer(firevect,{},'firevect'); 
  //Map.addLayer(firevect_ndvi,{},'fire_vect_NDVI');

// Make the histogram, set the options.
    //var histogram = ui.Chart.image.histogram(sub_nbr, region, 30);
 // print(histogram);
  