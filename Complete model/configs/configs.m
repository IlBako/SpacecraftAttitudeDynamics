astro_data = astronomicData;
sc_data = spaceCraftData;
orbit_data = orbitData;
in_cond = initialConditions(astro_data.n_eth);
pert_data = perturbationData(orbit_data.n_sc);
sensor_data = sensorData;
worldMag_data = load("data\WMD.mat");