A toy hydrologic model that uses vapor pressure, precipitation, lat/lon, elevation, incoming shortwave, temperature, and wind speed datasets to estiamte extraterrestrial and clear sky 
solar radiation, net longwave radiation, net radiation, PET, and ET (as a function of simulated PET and soil moisture). Finally, it employs 
fast/unsaturated/saturated reservoirs to simulate total streamflow. Details can be found in the docs section.

The main difference between this model (SVAT) and the Hydro_T model is - the Hydro_T model employs simplified Penman–Monteith (P-M) to estimate PET, while the SVAT model employs the full-scale P-M equation to estimate ET. Also, representation of the soil moisture is simple in the SVAT model - based on a combination of canopy storage, maximum canopy storage (const. for a land cover type) and canopy drainge (TH by SHuttleworth Fig 24.1 schematic). The SVAT model uses a combination of maximum cacopy storage, canopy storage, and canopy drainage to simulate soil moisture at every stage. However, the Hydro_T model employ a series of reservoirs to simulate soil moisture (starts with infiltration estimation using a fast flow resrvoir, initial storage and recharge calculation using an unsaturated zone reservoir, and finally the balance between infiltration, recharge, and maximum storage capacity to estimate current soil moisture).

__References:__ 
1. [Terrestrial Hydrometeorology by W. James Shuttleworth](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119951933)
2. Class notes


__Plots: Catchment 8__

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/ba172ad7-ab4c-4936-9538-681f4362ebf2)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/fa527826-5b17-4a1c-b2b1-19ca8fa5fbe6)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/38be03f3-67e1-491b-aa01-7bc7e406d185)











