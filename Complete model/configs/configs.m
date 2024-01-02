rng("default")

astro_data = astronomicData;
sc_data = spaceCraftData;
orbit_data = orbitData;
in_cond = initialConditions(orbit_data);
pert_data = perturbationData;
sensor_data = sensorData;
actuator_data = actuatorData;
worldMag_data = load("data\WMD.mat");