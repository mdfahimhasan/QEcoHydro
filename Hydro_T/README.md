A toy hydrologic model that uses vapor pressure, precipitation, lat/lon, elevation, incoming shortwave, temperature, and wind speed datasets to estiamte extraterrestrial and clear sky 
solar radiation, net longwave radiation, net radiation, PET, and ET (as a function of simulated PET and soil moisture). Finally, it employs 
fast/unsaturated/saturated reservoirs to simulate total streamflow. Details can be found in the docs section.

The main difference between this model (Hydro_T) and the SVAT model in this repository is - the Hydro_T model estimates PET as function of soil moisture (soil moisture availability ratio estimated from unsaturated zone storage and unsaturated zone maximus storage capacity), while the SVAT model employs the Penman–Monteith equation to estimate ET. 

__References:__ 
1. [Terrestrial Hydrometeorology by W. James Shuttleworth](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119951933)
2. Class notes


__Plots: Catchment 8__

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/ba172ad7-ab4c-4936-9538-681f4362ebf2)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/fa527826-5b17-4a1c-b2b1-19ca8fa5fbe6)

![image](https://github.com/mdfahimhasan/QEcoHydro/assets/77580408/38be03f3-67e1-491b-aa01-7bc7e406d185)











