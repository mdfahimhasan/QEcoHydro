Simple SVAT (Soil Vegetation Atmosphere Transfer Models) model generating ET and soil moisture outputs. Includes aerodynamic, surface resistance calculation model. Also, it has canopy water balance module to estimate ET (ET estimated as function of soil moisture; soil moisture availability ratio estimated from unsaturated zone storage and unsaturated zone maximus storage capacity) and simple soil moisture update module. The model has been prepared for the Quantitative Ecohydrology class (CIVE 625) as a course requirement. Currently, it can only simulate ET components and soil moisture for forest and grass type canipies. 

The main difference between this model (SVAT) and the Hydro_T model in this repository is - the Hydro_T model employs simplified Penman–Monteith (P-M) to estimate PET, while the SVAT model employs the full-scale P-M equation to estimate ET. Also, representation of the soil moisture is simple in the SVAT model - based on a combination of canopy storage, maximum canopy storage (const. for a land cover type) and canopy drainge (TH by SHuttleworth Fig 24.1 schematic). 


__Reference book:__ [Terrestrial Hydrometeorology by W. James Shuttleworth](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119951933)


__Acknowledgement:__ Thornton-Dunwoody, Alex (@ atdunwoody) for contributing in many functions used in different modules.


__Plots: for Grass canopy__

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/7b5b5cf8-4654-424f-8c00-ef76ea2a4325)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/96adb4e6-9864-473c-9c19-7bafb442dd0f)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/89f74520-209f-409c-80bc-cefca03cc271)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/867d9b51-db2d-4441-b008-591ed78f8ee0)








