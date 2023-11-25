"""
  ET CALCULATOR FOR LANDSAT 8 SATELLITE IMAGERY
  
  Updated on: 17 September 2022
  Written by HaKaN OGUZ
  Description: This python code calculates evapotranspiration (ET) from Landsat 8 satellite imagery 
"""

import arcpy
from sys import argv

def ETCalculatorfromLandsat8(band_2="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B2.TIF", band_3="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B3.TIF", band_4="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B4.TIF", band_5="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B5.TIF", band_6="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B6.TIF", band_7="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B7.TIF", band_10="C:\\ET_Calculator\\Landsat_8\\Input_Data\\LC08_L1TP_175034_20130913_20200912_02_T1_B10.TIF", SUN_ELEVATION=55.11, ESUN_DISTANCE=1.004, Upwelling_Path_Radiance=0.84, Atmospheric_Transmissivity=0.89, ET_Inst="C:\\ET_Calculator\\Landsat_8\\Output\\et_inst"):  # ET Calculator from Landsat 8

    # To allow overwriting outputs change overwriteOutput option to True.
    arcpy.env.overwriteOutput = False

    # Check out any necessary licenses.
    arcpy.CheckOutExtension("spatial")

    RADIANCE_MULT_BAND_10 = 0.000334
    RADIANCE_ADD_BAND_10 = 0.1
    Downwelling_Path_Radiance = 1.44
    Cold_Pixel_Temp_in_C = 11
    Wind_Speed_Height_m_ = 100
    Wind_Speed_m_s_ = 10
    RADIANCE_MULT_BAND_2 = 0.01275
    RADIANCE_ADD_BAND_2 = -63.78
    RADIANCE_MULT_BAND_3 = 0.01175
    RADIANCE_ADD_BAND_3 = -58.775
    RADIANCE_MULT_BAND_4 = 0.0099125
    RADIANCE_ADD_BAND_4 = -49.56261
    RADIANCE_MULT_BAND_5 = 0.006066
    RADIANCE_ADD_BAND_5 = -30.32984
    RADIANCE_MULT_BAND_6 = 0.001508
    RADIANCE_ADD_BAND_6 = -7.5427
    RADIANCE_MULT_BAND_7 = 0.000508
    RADIANCE_ADD_BAND_7 = -2.542
    Coefficient_b = -59
    Coefficient_a = 0.204

    # Process: Calculate Radiance for band 10 (Raster Calculator) 
    band10_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band10_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_10, band_10, RADIANCE_ADD_BAND_10], output_raster=band10_radiance)

    # Process: Calculate Black Body Temp (Raster Calculator) 
    Tb_temp = "C:\\ET_Calculator\\Landsat_8\\Output\\tb_temp"
    arcpy.gp.RasterCalculator_sa(expression=[band10_radiance], output_raster=Tb_temp)

    # Process: Calculate Value for Friction Velocity (U*) (Calculate Value) 
    Friction_Velocity_U_ = (0.41 * float(Wind_Speed_m_s_)) / (math.log(float(Wind_Speed_Height_m_) / (0.12 * 0.3)))[0]

    # Process: Calculate Radiance for band 2 (Raster Calculator) 
    band2_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band2_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_2, band_2, RADIANCE_ADD_BAND_2], output_raster=band2_radiance)

    # Process: Calculate Plenetary Reflectance for band 2 (Raster Calculator) 
    band2_reflectance = "C:\\ET_Calculator\\Landsat_8\\Output\\band2_ref"
    arcpy.gp.RasterCalculator_sa(expression=[band2_radiance, SUN_ELEVATION, ESUN_DISTANCE], output_raster=band2_reflectance)

    # Process: Calculate Radiance for band 3 (Raster Calculator) 
    band3_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band3_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_3, band_3, RADIANCE_ADD_BAND_3], output_raster=band3_radiance)

    # Process: Calculate Plenetary Reflectance for band 3 (Raster Calculator) 
    band3_reflectance = "C:\\ET_Calculator\\Landsat_8\\Output\\band3_ref"
    arcpy.gp.RasterCalculator_sa(expression=[band3_radiance, SUN_ELEVATION, ESUN_DISTANCE], output_raster=band3_reflectance)

    # Process: Calculate Radiance for band 4 (Raster Calculator) 
    band4_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band4_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_4, band_4, RADIANCE_ADD_BAND_4], output_raster=band4_radiance)

    # Process: Calculate Plenetary Reflectance for band 4 (Raster Calculator) 
    band4_reflectance = "C:\\ET_Calculator\\Landsat_8\\Output\\band4_ref"
    arcpy.gp.RasterCalculator_sa(expression=[band4_radiance, SUN_ELEVATION, ESUN_DISTANCE], output_raster=band4_reflectance)

    # Process: Calculate Radiance for band 5 (Raster Calculator) 
    band5_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band5_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_5, band_5, RADIANCE_ADD_BAND_5], output_raster=band5_radiance)

    # Process: Calculate Plenetary Reflectance for band 5 (Raster Calculator) 
    band5_reflectance = "C:\\ET_Calculator\\Landsat_8\\Output\\band5_ref"
    arcpy.gp.RasterCalculator_sa(expression=[band5_radiance, SUN_ELEVATION, ESUN_DISTANCE], output_raster=band5_reflectance)

    # Process: Calculate Radiance for band 6 (Raster Calculator) 
    band6_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band6_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_6, band_6, RADIANCE_ADD_BAND_6], output_raster=band6_radiance)

    # Process: Calculate Plenetary Reflectance for band 6 (Raster Calculator) 
    band6_reflectance = "C:\\ET_Calculator\\Landsat_8\\Output\\band6_ref"
    arcpy.gp.RasterCalculator_sa(expression=[band6_radiance, SUN_ELEVATION, ESUN_DISTANCE], output_raster=band6_reflectance)

    # Process: Calculate Radiance for band 7 (Raster Calculator) 
    band7_radiance = "C:\\ET_Calculator\\Landsat_8\\Output\\band7_rad"
    arcpy.gp.RasterCalculator_sa(expression=[RADIANCE_MULT_BAND_7, band_7, RADIANCE_ADD_BAND_7], output_raster=band7_radiance)

    # Process: Calculate Plenetary Reflectance for band 7 (Raster Calculator) 
    band7_reflectance = "C:\\ET_Calculator\\Landsat_8\\Output\\band7_ref"
    arcpy.gp.RasterCalculator_sa(expression=[band7_radiance, SUN_ELEVATION, ESUN_DISTANCE], output_raster=band7_reflectance)

    # Process: Calculate Albedo TOA (Raster Calculator) 
    Albedo_toa = "C:\\ET_Calculator\\Landsat_8\\Output\\albedo_toa"
    arcpy.gp.RasterCalculator_sa(expression=[band2_reflectance, band3_reflectance, band4_reflectance, band5_reflectance, band6_reflectance, band7_reflectance], output_raster=Albedo_toa)

    # Process: Calculate Surface Albedo (Raster Calculator) 
    Surface_Albedo = "C:\\ET_Calculator\\Landsat_8\\Output\\surf_albedo"
    arcpy.gp.RasterCalculator_sa(expression=[Albedo_toa, Atmospheric_Transmissivity, Atmospheric_Transmissivity], output_raster=Surface_Albedo)

    # Process: Calculate Incoming Shortwave Radiation (Calculate Value) 
    RS_Down = float(Atmospheric_Transmissivity) * (1367 / (float(ESUN_DISTANCE) * float(ESUN_DISTANCE))) * math.cos((90 - float(SUN_ELEVATION)) * math.pi / 180)[0]

    # Process: Calculate Incoming Longwave Radiation (Calculate Value) 
    RL_Down = (0.85 * (0 - math.log(float(Atmospheric_Transmissivity))) ** 0.09) * (5.67 / 100000000) * ((float(Cold_Pixel_Temp_in_C) + 273.15) ** 4)[0]

    # Process: Calculate NDVI (Raster Calculator) 
    NDVI = "C:\\ET_Calculator\\Landsat_8\\Output\\ndvi"
    arcpy.gp.RasterCalculator_sa(expression=[band5_reflectance, band4_reflectance, band5_reflectance, band4_reflectance], output_raster=NDVI)

    # Process: Calculate NDVI Corrected (Raster Calculator) 
    NDVI_Corrected = "C:\\ET_Calculator\\Landsat_8\\Output\\ndvi_correct"
    arcpy.gp.RasterCalculator_sa(expression=[NDVI, NDVI, NDVI, NDVI, NDVI], output_raster=NDVI_Corrected)

    # Process: Calculate LAI (Raster Calculator) 
    LAI = "C:\\ET_Calculator\\Landsat_8\\Output\\lai"
    arcpy.gp.RasterCalculator_sa(expression=[NDVI_Corrected], output_raster=LAI)

    # Process: Calculate Emissivity (Raster Calculator) 
    Emissivity = "C:\\ET_Calculator\\Landsat_8\\Output\\emissivity"
    arcpy.gp.RasterCalculator_sa(expression=[NDVI_Corrected, LAI, NDVI_Corrected, LAI, LAI, NDVI_Corrected, Surface_Albedo, NDVI_Corrected, Surface_Albedo], output_raster=Emissivity)

    # Process: Calculate Rc (Raster Calculator) 
    Rc = "C:\\ET_Calculator\\Landsat_8\\Output\\rc"
    arcpy.gp.RasterCalculator_sa(expression=[band10_radiance, Upwelling_Path_Radiance, Atmospheric_Transmissivity], output_raster=Rc)

    # Process: Calculate Surface Temp (Raster Calculator) 
    Surface_Temperature = "C:\\ET_Calculator\\Landsat_8\\Output\\surf_temp"
    arcpy.gp.RasterCalculator_sa(expression=[Emissivity, Rc], output_raster=Surface_Temperature)

    # Process: Calculate Outgoing Longwave Radiation (Raster Calculator) 
    RL_up = "C:\\ET_Calculator\\Landsat_8\\Output\\rl_up"
    arcpy.gp.RasterCalculator_sa(expression=[Emissivity, Surface_Temperature], output_raster=RL_up)

    # Process: Calculate Net Surface Radiation (Rn) (Raster Calculator) 
    Rn = "C:\\ET_Calculator\\Landsat_8\\Output\\rn"
    arcpy.gp.RasterCalculator_sa(expression=[Surface_Albedo, RS_Down, RL_Down, RL_up, Emissivity, RL_Down], output_raster=Rn)

    # Process: Calculate G/Rn (Raster Calculator) 
    G_Rn = "C:\\ET_Calculator\\Landsat_8\\Output\\g_rn"
    arcpy.gp.RasterCalculator_sa(expression=[Surface_Temperature, Surface_Albedo, Surface_Albedo, Surface_Albedo, Surface_Albedo, NDVI_Corrected], output_raster=G_Rn)

    # Process: Calculate G/Rn Corrected (Raster Calculator) 
    G_Rn_Corrected = "C:\\ET_Calculator\\G_Rn_correct"
    arcpy.gp.RasterCalculator_sa(expression=[NDVI_Corrected, Surface_Temperature, Surface_Albedo, G_Rn], output_raster=G_Rn_Corrected)

    # Process: Calculate G (Raster Calculator) 
    G = "C:\\ET_Calculator\\Landsat_8\\Output\\g"
    arcpy.gp.RasterCalculator_sa(expression=[Rn, G_Rn_Corrected], output_raster=G)

    # Process: Calculate dT (Raster Calculator) 
    dT = "C:\\ET_Calculator\\Landsat_8\\Output\\dt"
    arcpy.gp.RasterCalculator_sa(expression=[Coefficient_b, Coefficient_a, Surface_Temperature], output_raster=dT)

    # Process: Calculate U_200 (Calculate Value) 
    u_200 = float(Friction_Velocity_U_) * ((math.log(float(Wind_Speed_Height_m_) / (0.12 * 0.3))) / 0.41)[0]

    # Process: Calculate Friction Velocity U*_new (Raster Calculator) 
    U_new = "C:\\ET_Calculator\\Landsat_8\\Output\\u_new"
    arcpy.gp.RasterCalculator_sa(expression=[u_200, LAI], output_raster=U_new)

    # Process: Calculate Rah (Raster Calculator) 
    Rah = "C:\\ET_Calculator\\Landsat_8\\Output\\rah"
    arcpy.gp.RasterCalculator_sa(expression=[U_new], output_raster=Rah)

    # Process: Calculate H (Raster Calculator) 
    H = "C:\\ET_Calculator\\Landsat_8\\Output\\h"
    arcpy.gp.RasterCalculator_sa(expression=[dT, Rah], output_raster=H)

    # Process: Calculate Latent Heat Flux (Latent ET) (Raster Calculator) 
    Latent_ET = "C:\\ET_Calculator\\Landsat_8\\Output\\latent_ET"
    arcpy.gp.RasterCalculator_sa(expression=[Rn, G, H], output_raster=Latent_ET)

    # Process: Calculate ET Inst (Raster Calculator) 
    arcpy.gp.RasterCalculator_sa(expression=[Latent_ET], output_raster=ET_Inst)

if __name__ == '__main__':
    # Global Environment settings
    with arcpy.EnvManager(scratchWorkspace=r"C:\Temp\MyProject1\MyProject1.gdb", workspace=r"C:\Temp\MyProject1\MyProject1.gdb"):
        ETCalculatorfromLandsat8(*argv[1:])
